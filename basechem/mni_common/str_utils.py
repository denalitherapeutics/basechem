import re


def capitalize_first(string):
    """
    Capitalize the first letter of a string. This is used to capitalize the first letter of
    verbose names in the `get_choices_for_field` function
    """
    return string[0].upper() + string[1:]


def natsort(string):
    """
    Calculates the natural sort key for a string. Natural sorting is helpful when sorting strings
    that include numbers, ex: "a1", "A2", ..., "a10", "A11", since traditional sort will put numbers
    like "10" before "2". Use this function like this sorted(list, key=natsort)
    :param string: a string to calculate the natsort key for
    """
    return [int(c) if c.isdigit() else c.lower() for c in re.split("(\d+)", string)]


def get_A1_column_name_from_index(col_idx):
    """
    Convert a column index to a column name using the A1 reference style
    Examples:
        1 -> "A"
        26 -> "Z"
        27 -> "AA"
        52 -> "AZ"
    :param col_idx: an integer, the column index to convert
    :return: a string, the A1 style column name
    """
    if not isinstance(col_idx, int) or col_idx < 1:
        raise ValueError(f"Invalid column index: {col_idx}")
    letters = []
    while col_idx > 0:
        col_idx, remainder = divmod(col_idx, 26)
        if remainder == 0:
            remainder = 26
            col_idx -= 1
        letters.append(chr(remainder + 64))
    return "".join(reversed(letters))


def get_A1_column_index_from_name(col_name):
    """
    Convert a column string into a column index using the A1 reference style
    Examples:
        "A" -> 1
        "Z" -> 26
        "AA" -> 27
        "AZ" -> 52
    :param col_name: a string, the Excel-style column name to convert
    :return: an integer, the column index
    """
    reversed_letters = col_name.strip()[::-1]
    col_index = 0
    factor = 1
    for letter in reversed_letters:
        letter = letter.upper()
        if not "A" <= letter <= "Z":
            raise ValueError(f"Invalid column letter {letter}")
        col_index += (ord(letter) - 64) * factor
        factor *= 26
    return col_index
