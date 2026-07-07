"""Domain exceptions for the VCF reader."""


class VCFError(Exception):
    """Base exception for VCF reader errors."""


class VCFValidationError(VCFError):
    """Raised when a VCF file fails structural validation."""


class VCFParseError(VCFError):
    """Raised when a VCF line cannot be parsed."""


class VCFExportError(VCFError):
    """Raised when an export operation fails."""


class VCFReaderStateError(VCFError):
    """Raised when reader methods are called in an invalid state."""
