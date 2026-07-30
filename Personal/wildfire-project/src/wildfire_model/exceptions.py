"""Domain-specific errors with explicit failure semantics."""


class WildfireModelError(Exception):
    """Base class for all model errors."""


class InvalidConfigurationError(WildfireModelError, ValueError):
    """Raised when a physical or numerical configuration is invalid."""


class NonPhysicalStateError(WildfireModelError, ValueError):
    """Raised when temperature or biomass leaves its physical domain."""


class IntegrationConvergenceError(WildfireModelError, RuntimeError):
    """Raised when an implicit integration stage does not converge."""
