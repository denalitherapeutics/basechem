from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D


def moltext_to_svg(moltext, size_x=450, size_y=150, transparent=True):
    """
    Given a string of moltext, generate an SVG string with the given dimensions depicting this molecule
    :param moltext: a string, valid moltext that should be drawn as an SVG
    :param size_x: an integer, the width of the returned svg (in pixels)
    :param size_y: an integer, the height of the returned svg (in pixels)
    :param transparent: a boolean, should the background be transparent? If False, the background is white
    :returns: a string with SVG data
    """
    mol = Chem.MolFromMolBlock(moltext)
    if not mol:
        return ""
    try:
        mol.GetAtomWithIdx(0).GetExplicitValence()
    except RuntimeError:
        mol.UpdatePropertyCache(False)
    try:
        mc_mol = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=True)
    except ValueError:  # <- can happen on a kekulization failure
        mc_mol = rdMolDraw2D.PrepareMolForDrawing(mol, kekulize=False)
    # Initialize the drawer with the correct canvas size
    drawer = rdMolDraw2D.MolDraw2DSVG(size_x, size_y)
    # Update settings like bond width and font size
    options = drawer.drawOptions()
    options.bondLineWidth = 1
    options.maxFontSize = 18
    options.fixedBondLength = 50
    if transparent:
        options.clearBackground = False
    drawer.SetDrawOptions(options)
    # Draw the molecule
    drawer.DrawMolecule(mc_mol)
    drawer.FinishDrawing()
    svg = drawer.GetDrawingText()
    return svg
