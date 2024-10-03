import re

from django import template
from django.utils.safestring import mark_safe

from basechem.main.models.compound_models import Compound, Series
from basechem.main.utils import strip_series_from_conf_id

register = template.Library()


@register.filter(name="as_svg")
def as_svg(compound, dimensions):
    """
    Generate an SVG image with a white background for the given compound
    :param compound: a Compound object
    :param dimensions: a string of the form "x_dim, y_dim"
    :returns: an svg with the 2d structure of the compound
    """
    [x, y] = dimensions.split(",")
    return compound.get_inline_svg(int(x), int(y))


@register.filter(name="as_transparent_svg")
def as_transparent_svg(compound, dimensions):
    """
    Generate an SVG image with a transparent background for the given compound
    :param compound: a Compound object
    :param dimensions: a string of the form "x_dim, y_dim"
    :returns: an svg with the 2d structure of the compound
    """
    [x, y] = dimensions.split(",")
    return compound.get_inline_svg(int(x), int(y), transparent=True)


@register.filter(name="get_collection_urls")
def get_collection_urls(compound):
    """
    This function is called from the `project_home.html` template used for the `ProjectView` to
    display a list of collections that use a particular compound.
    :param compound: a Compound object
    :return: a list of url paths for each collection that contains the compound
    """
    # Because this function is called once for each compound in the project (a lot of calls),
    # it has been specifically written to not make an additional DB call. `compoundoccurrence_set`
    # and `collection_set` are cached in the queryset that `compound` is coming from
    coll_ids = set()
    for co in compound.compoundoccurrence_set.all():
        coll_ids.update(coll.pk for coll in co.collection_set.all())
    urls = []
    for cid in sorted(coll_ids):
        url = f"<a href=/align/{cid}>{cid}</a>"
        urls.append(url)

    return mark_safe(", ".join(urls))


@register.filter(name="get_most_recent_date")
def get_most_recent_date(compound):
    """
    :param compound: a Compound object
    :return: the most recent generated date for all COs for this Compound
    """
    try:
        # Most recent CO. Converting to list because doing `.last()` is a new query
        # that won't used the cached `compoundoccurrence_set.all()`
        co = list(compound.compoundoccurrence_set.all())[-1]
        return co.generated.strftime("%b. %d, %Y")
    except:
        return None


@register.filter(name="get_user_ids")
def get_user_ids(compound):
    """
    This function is called from the `project_home.html` template used for the `ProjectView` to
    display a list of users who have a collection with `compound`
    :param compound: a Compound object
    :return: a string, comma separated list of users who have collections with the given compound
    """
    # Because this function is called once for each compound in the project (a lot of calls),
    # it has been specifically written to not make an additional DB call. `compoundoccurrence_set`
    # is cached in the queryset that `compound` is coming from
    cos = compound.compoundoccurrence_set.all()
    users = [co.owner.last_name for co in cos]
    return ", ".join(set(users))


@register.filter(name="get_conformers")
def get_conformers(data, co_id):
    """
    :param data: a dictionary of moltext returned by the ligand alignment backend
    :param co_id: an int, the ID of the compound occurrence
    :return: a list of tuples of the form (conformer ID, moltext)
    """
    str_co_id = f"co-{co_id}"
    if data:
        return data["compounds"].get(str_co_id, {}).items()
    else:
        return []


@register.filter(name="get_display_name")
def get_display_name(reference):
    """
    :param reference: a string 'c-pk' or 's-pk' referring to the reference compound
    :return: a display name for the reference (ex. 'Compound: DN00001')
    """
    reference_id = int(reference[2:])
    if "c-" in reference:
        compound = Compound.objects.get(id=reference_id)
        return f"Compound: {compound.name}"
    elif "s-" in reference:
        series = Series.objects.get(id=reference_id)
        return f"Series: {series.name}"


@register.filter(name="color_filter")
def color_filter(value, prop):
    """
    :param value: the value to determine what color it should be
    :param prop: the property to color
    :return: the color the value should be colored
    """
    # Prediction Colors
    if prop == "HLM Prediction":
        value = float(value)
        if value <= 6.2:
            return "green"
        elif 6.2 < value <= 15.0:
            return "orange"
        elif value > 15.0:
            return "red"

    if prop == "RLM Prediction":
        value = float(value)
        if value <= 17:
            return "green"
        elif 17 < value <= 39.0:
            return "orange"
        elif value > 39.0:
            return "red"

    if prop == "TPSA":
        value = float(value)
        if value < 100:
            return "green"
        elif 100 <= value < 120:
            return "orange"
        elif value >= 120:
            return "red"

    if prop == "LogD Prediction":
        value = float(value)
        if value > 4 or value < 1:
            return "red"
        elif value > 3 or value < 2:
            return "orange"
        else:
            return "green"

    if prop == "assay_result" and value:
        value = float(value)
        if value <= 1:
            return "green"
        elif value > 10:
            return "red"
        else:
            return "orange"


@register.filter(name="get_delta_energy")
def get_delta_energy(data, co_id):
    """
    :param data: a dictionary
    :param co_id: a pk to get data for
    :return: the value of this key in the dictionary (or an empty list)
    """
    data = data.get(f"co-{co_id}", {})
    if data.get("delta_energy"):
        return "%.3f" % float(data["delta_energy"])
    return ""


@register.filter(name="get_propcalc_val")
def get_propcalc_val(prop_dict, prop):
    """
    Given a dictionary of properties and a propcalc prop, return the value as expected by the propcalc html template
    :param prop_dict: a dictionary of properties and their values as returned by `mol_w_properties` (lowercase, w/ underscores)
    :param prop: the user-facing property name (uppercase, no underscores)
    """
    prop = prop.lower().replace(" ", "_")
    value = prop_dict.get(prop, "")

    if "probabilities" in prop:
        prefix = re.match("(r|h)lm", prop).group(0)
        # Add the out of domain flag to the value
        ood_flag = eval(prop_dict.get(f"{prefix}_ood"))
        value = [value, ood_flag]

    if "prediction" in prop:
        prefix = prop
        interp_img = None

        match = re.match("(r|h)lm", prop)
        if match:
            prefix = re.match("(r|h)lm", prop).group(0)

        # Show measured value as 'value' if exists
        measured = prop_dict.get(f"{prefix}_measured")
        is_measured = False
        if measured not in ["", None]:
            is_measured = True
            value = measured

        # Add the interpretation image path to the value
        interp_img = prop_dict.get(f"{prefix}_interpretation")
        if interp_img:
            value = [value, is_measured, interp_img]
        else:
            value = [value, is_measured]
    return value


@register.filter(name="pred_type")
def pred_type(value):
    """
    Determine if the value is measured or predicted
    :param value: a list of [str(pred), is_measured...)]
    :return: "measured" if there is a measured value
    """
    if type(value) != list:
        return ""
    if len(value) >= 2:
        if value[1]:  # Expected location for `measured` bool
            return "measured"
    return ""


@register.filter(name="conf_id_display")
def conf_id_display(conf_id):
    """
    :param conf_id: a string, the full conf_id from an analysis (ex: "c120-co590-1-1-s-2")
    """
    return strip_series_from_conf_id(conf_id)


@register.filter(name="compared_color")
def compared_color(first_val, second_val):
    """
    :param first_val: first value to compare
    :param second_val: second value to compare
    :return: a string, "green" if the first value is less than the second value, "red" if not.
    """
    if first_val and second_val:
        if float(first_val) < float(second_val):
            return "green"
        else:
            return "red"
