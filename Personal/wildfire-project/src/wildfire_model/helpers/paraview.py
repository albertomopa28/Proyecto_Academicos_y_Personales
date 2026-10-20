"""XDMF/HDF5 time series for ParaView, with cell data on a terrain surface."""

from pathlib import Path
from xml.etree import ElementTree as ET

import h5py
import numpy as np
from numpy.typing import NDArray

from wildfire_model.exceptions import InvalidConfigurationError
from wildfire_model.inputs.fields import FieldInput, evaluate_field
from wildfire_model.mesh.mesh_1d import UniformMesh1D
from wildfire_model.mesh.mesh_2d import UniformMesh2D
from wildfire_model.models.state import validate_temperature_biomass
from wildfire_model.numerics.solver_1d import SimulationResult
from wildfire_model.numerics.solver_2d import SimulationResult2D


def _cell_elevation_to_vertices(values: NDArray[np.float64]) -> NDArray[np.float64]:
    """Interpolate cell heights to vertices, with linear edge extrapolation."""

    result = values
    for axis in (1, 0):
        moved = np.moveaxis(result, axis, -1)
        vertices = np.empty((*moved.shape[:-1], moved.shape[-1] + 1), dtype=np.float64)
        if moved.shape[-1] == 1:
            vertices[..., :] = moved
        else:
            vertices[..., 1:-1] = 0.5 * (moved[..., :-1] + moved[..., 1:])
            vertices[..., 0] = 1.5 * moved[..., 0] - 0.5 * moved[..., 1]
            vertices[..., -1] = 1.5 * moved[..., -1] - 0.5 * moved[..., -2]
        result = np.moveaxis(vertices, -1, axis)
    return result


def _terrain_points(mesh: UniformMesh2D, elevation: NDArray[np.float64]) -> NDArray[np.float64]:
    x, y = np.meshgrid(mesh.x_faces, mesh.y_faces, indexing="xy")
    z = _cell_elevation_to_vertices(elevation)
    return np.column_stack((x.ravel(), y.ravel(), z.ravel()))


def _quadrilaterals(mesh: UniformMesh2D) -> NDArray[np.int64]:
    nodes = np.arange((mesh.ny + 1) * (mesh.nx + 1), dtype=np.int64).reshape(
        mesh.ny + 1, mesh.nx + 1
    )
    return np.column_stack(
        (
            nodes[:-1, :-1].ravel(),
            nodes[:-1, 1:].ravel(),
            nodes[1:, 1:].ravel(),
            nodes[1:, :-1].ravel(),
        )
    )


def _data_item(
    parent: ET.Element,
    shape: tuple[int, ...],
    reference: str,
    *,
    integer: bool = False,
) -> None:
    item = ET.SubElement(
        parent,
        "DataItem",
        Dimensions=" ".join(str(value) for value in shape),
        DataType="Int" if integer else "Float",
        Precision="8",
        Format="HDF",
    )
    item.text = reference


def write_paraview_series(
    mesh: UniformMesh1D | UniformMesh2D,
    result: SimulationResult | SimulationResult2D,
    *,
    output_directory: str | Path = "downloads/simulation",
    elevation: FieldInput = 0.0,
) -> Path:
    """Write ``series.xdmf`` and ``series.h5`` under a downloads directory.

    Fields retain their finite-volume cell association and SI units. A 2D
    mesh is drawn at its supplied elevation; a 1D mesh is drawn along x.
    Terrain vertices are interpolated from the cell field used by the solver.
    """

    times = np.asarray(result.times, dtype=np.float64)
    if (
        times.ndim != 1
        or times.size == 0
        or not np.all(np.isfinite(times))
        or times[0] < 0.0
        or np.any(np.diff(times) <= 0.0)
    ):
        raise InvalidConfigurationError("Export times must be finite, non-negative and increasing.")
    shape = mesh.shape if isinstance(mesh, UniformMesh2D) else (mesh.n_cells,)
    expected = (times.size, *shape)
    temperature, biomass = validate_temperature_biomass(
        result.temperature, result.biomass, context="ParaView export"
    )
    if temperature.shape != expected:
        raise InvalidConfigurationError(f"Export fields must have shape {expected}.")
    directory = Path(output_directory).expanduser().resolve()
    if "downloads" not in directory.parts:
        raise InvalidConfigurationError("ParaView results must be written under downloads/.")

    if isinstance(mesh, UniformMesh2D):
        cells = _quadrilaterals(mesh)
        topology_type = "Quadrilateral"
    else:
        cells = np.column_stack(
            (
                np.arange(mesh.n_cells, dtype=np.int64),
                np.arange(1, mesh.n_cells + 1, dtype=np.int64),
            )
        )
        topology_type = "Polyline"
    directory.mkdir(parents=True, exist_ok=True)
    xdmf_path = directory / "series.xdmf"
    hdf_path = directory / "series.h5"
    root = ET.Element("Xdmf", Version="3.0")
    domain = ET.SubElement(root, "Domain")
    collection = ET.SubElement(
        domain, "Grid", Name="Wildfire", GridType="Collection", CollectionType="Temporal"
    )
    previous_elevation: NDArray[np.float64] | None = None
    points_path = ""
    with h5py.File(hdf_path, "w") as data:
        data.attrs["model"] = "Vertically averaged wildfire; terrain is a visualization surface"
        data.attrs["cell_shape"] = shape
        data.create_dataset("mesh/cells", data=cells)
        data.create_dataset("times_s", data=times)
        for index, time in enumerate(times):
            if isinstance(mesh, UniformMesh2D):
                heights = evaluate_field(
                    elevation, mesh, name="terrain elevation", time=float(time)
                )
                if previous_elevation is None or not np.array_equal(heights, previous_elevation):
                    points = _terrain_points(mesh, heights)
                    points_path = f"mesh/points_{index:06d}"
                    data.create_dataset(points_path, data=points, compression="gzip")
                    previous_elevation = heights
            else:
                heights = np.zeros(mesh.n_cells, dtype=np.float64)
                points = np.column_stack((mesh.faces, np.zeros((mesh.n_cells + 1, 2))))
                points_path = "mesh/points"
                if index == 0:
                    data.create_dataset(points_path, data=points)

            grid = ET.SubElement(collection, "Grid", Name=f"step_{index:06d}", GridType="Uniform")
            ET.SubElement(grid, "Time", Value=format(float(time), ".17g"))
            topology = ET.SubElement(
                grid, "Topology", TopologyType=topology_type, NumberOfElements=str(cells.shape[0])
            )
            if topology_type == "Polyline":
                topology.set("NodesPerElement", "2")
            _data_item(topology, cells.shape, f"{hdf_path.name}:/mesh/cells", integer=True)
            geometry = ET.SubElement(grid, "Geometry", GeometryType="XYZ")
            _data_item(geometry, points.shape, f"{hdf_path.name}:/{points_path}")
            fields = {
                "temperature_K": temperature[index],
                "remaining_fuel": biomass[index],
                "burned_fraction": 1.0 - biomass[index],
                "elevation_m": heights,
            }
            for name, values in fields.items():
                dataset_path = f"steps/{index:06d}/{name}"
                data.create_dataset(dataset_path, data=values.ravel(), compression="gzip")
                attribute = ET.SubElement(
                    grid, "Attribute", Name=name, AttributeType="Scalar", Center="Cell"
                )
                _data_item(attribute, (cells.shape[0],), f"{hdf_path.name}:/{dataset_path}")
    ET.indent(root)
    ET.ElementTree(root).write(xdmf_path, encoding="utf-8", xml_declaration=True)
    return xdmf_path
