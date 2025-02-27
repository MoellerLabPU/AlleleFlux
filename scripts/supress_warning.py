# utilities.py
import logging
import warnings
from collections import defaultdict

# A dictionary to count occurrences for each warning type.
_warning_counter = defaultdict(int)

# List of warning message substrings to check.
SUPPRESS_MESSAGES = [
    "Sample size too small for normal approximation",
    "Precision loss occurred in moment calculation due to catastrophic cancellation.",
]


def custom_showwarning(message, category, filename, lineno, file=None, line=None):
    """
    Custom warning handler to control the display of specific warning messages.

    This function intercepts warning messages and checks if they contain the
    phrase "Sample size too small for normal approximation". If the warning
    message matches this phrase, it counts the occurrences and only displays
    the warning up to two times. After the second occurrence, it suppresses
    further warnings of this type and logs a message indicating that further
    warnings will be suppressed.

    Since we are doing multiprocessing, each process will track its own _warning_counter,
    meaning we might see up to 2 * num_processes warnings.

    For other warning messages, it displays them normally.

    Parameters:
    - message (Warning): The warning message being processed.
    - category (Warning): The category of the warning.
    - filename (str): The name of the file where the warning occurred.
    - lineno (int): The line number where the warning occurred.
    - file (file object, optional): The file to which the warning should be written.
    - line (str, optional): The source code line that triggered the warning.

    Returns:
    None
    """
    msg_str = str(message)
    # Check if the message contains any of the target phrases.
    for phrase in SUPPRESS_MESSAGES:
        if phrase in msg_str:
            key = (phrase, category)
            _warning_counter[key] += 1
            if _warning_counter[key] <= 2:
                warnings._showwarning_orig(
                    message, category, filename, lineno, file, line
                )
                if _warning_counter[key] == 2:
                    logging.warning(
                        f"Further warnings matching '{phrase}' will be suppressed.\n"
                    )
            # Return early so that warnings matching our criteria aren't processed further.
            return
    # For all other warnings, show them normally.
    warnings._showwarning_orig(message, category, filename, lineno, file, line)


# Backup the original warning display function.
warnings._showwarning_orig = warnings.showwarning
# Override it with our custom function.
warnings.showwarning = custom_showwarning
