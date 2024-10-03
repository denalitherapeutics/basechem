import re


def strip_series_from_conf_id(conf_id):
    """
    When possible, remove the series component from a conformer ID
    :param conf_id: a string, a conformer/pose ID returned from a backend analysis
    :returns: a string, the conformer ID without the series
    """
    matches = re.match("(^.*)(-s-\d+)$", conf_id)
    if matches:
        return matches[1]
    else:
        return conf_id
