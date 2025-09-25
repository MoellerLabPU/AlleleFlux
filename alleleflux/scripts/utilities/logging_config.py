"""Centralized logging configuration utilities for AlleleFlux.

This module exposes a single function, `setup_logging`, that configures a
reasonably safe and consistent application-wide logging setup using
`logging.config.dictConfig`.

Behavior highlights:
- Honors the `ALLELEFLUX_LOG_LEVEL` environment variable when present.
- Accepts either numeric levels (e.g., "20") or level names (e.g., "INFO").
- Emits a concise warning to stderr and falls back to INFO on invalid values.
- Supports colored output for different log levels when terminal supports it.

The configuration attaches a single `StreamHandler` to stderr with a consistent
format that includes timestamp, level, logger name, and message, with optional
colors for better readability.
"""

import logging
import logging.config
import os
import sys
from typing import Union


class ColoredFormatter(logging.Formatter):
    """Custom formatter that adds ANSI color codes to log levels."""

    # ANSI color codes
    COLORS = {
        "DEBUG": "\033[36m",  # Cyan
        "INFO": "\033[37m",  # White (default)
        "WARNING": "\033[33m",  # Yellow
        "ERROR": "\033[31m",  # Red
        "CRITICAL": "\033[35m",  # Magenta
    }
    RESET = "\033[0m"  # Reset color

    def __init__(self, *args, use_color=None, **kwargs):
        super().__init__(*args, **kwargs)
        # Auto-detect color support if not explicitly set
        if use_color is None:
            # Check if we're outputting to a terminal that supports colors
            self.use_color = (
                hasattr(sys.stderr, "isatty")
                and sys.stderr.isatty()
                and os.getenv("TERM") != "dumb"
                and os.getenv("NO_COLOR") is None
            )
        else:
            self.use_color = use_color

    def format(self, record):
        # Get the original formatted message
        message = super().format(record)

        if self.use_color:
            # Add color to the level name in the message
            level_name = record.levelname
            if level_name in self.COLORS:
                # Replace the level name with colored version
                colored_level = f"{self.COLORS[level_name]}{level_name}{self.RESET}"
                message = message.replace(f"[{level_name}]", f"[{colored_level}]")

        return message


def setup_logging(
    level: Union[int, str] = logging.INFO, use_color: bool = None
) -> None:
    """Configure root logger via dictConfig with env-based level override and optional colors.

    Parameters
    ----------
    level : int | str, optional
        Default logging level to use when the `ALLELEFLUX_LOG_LEVEL` environment
        variable is not set. Can be an integer (e.g., ``logging.INFO`` or ``20``)
        or a case-insensitive level name such as ``"INFO"`` or ``"debug"``.
    use_color : bool, optional
        Whether to use colored output. If None (default), auto-detects based on
        terminal capabilities and environment variables (NO_COLOR, TERM).

    Notes
    -----
    - If `ALLELEFLUX_LOG_LEVEL` is set, its value takes precedence.
    - Invalid values result in a warning to stderr and a fallback to INFO.
    - Colors are automatically disabled if output is not a TTY or if NO_COLOR is set.
    - This function configures the root logger with one console handler that
      writes to stderr. Subsequent calls re-apply the same configuration; since
      `dictConfig` replaces the root configuration, duplicate handlers are not
      accumulated.
    """
    # Read level from environment if provided, falling back to the default `level`.
    raw = os.getenv("ALLELEFLUX_LOG_LEVEL", str(level)).strip()

    # Coerce level robustly:
    # - Numeric strings like "20" are converted to int
    # - Names like "info"/"WARNING" are mapped to logging constants
    # - Anything else triggers a warning and falls back to INFO
    try:
        eff_level = int(raw)
    except ValueError:
        eff_level = getattr(logging, raw.upper(), None)
        if not isinstance(eff_level, int):
            print(f"Warning: Invalid log level '{raw}', using INFO.", file=sys.stderr)
            eff_level = logging.INFO

    # Centralized logging configuration using dictConfig. This ensures a
    # predictable formatter and avoids accidental duplication of handlers that
    # can happen with successive `basicConfig` calls.
    config = {
        "version": 1,
        # Keep existing non-root loggers enabled so library loggers still work
        # unless explicitly disabled by the application.
        "disable_existing_loggers": False,
        "formatters": {
            "colored": {
                "()": ColoredFormatter,
                "format": "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
                "datefmt": "%Y-%m-%d %H:%M:%S",
                "use_color": use_color,
            },
            "plain": {
                "format": "%(asctime)s [%(levelname)s] %(name)s: %(message)s",
                "datefmt": "%Y-%m-%d %H:%M:%S",
            },
        },
        "handlers": {
            "console": {
                "class": "logging.StreamHandler",
                "formatter": "colored",
                # Send logs to stderr so stdout remains usable for data output
                "stream": "ext://sys.stderr",
                "level": eff_level,
            }
        },
        "root": {
            "level": eff_level,
            "handlers": ["console"],
        },
    }
    # Apply the configuration atomically.
    logging.config.dictConfig(config)
    logging.getLogger(__name__).info(
        f"Logging configured with level: {logging.getLevelName(eff_level)}"
    )
