/**
* Because we can't interact with ketcher until the iframe has loaded, this function
* waits for the iframe to load and then calls ketcherLoaded.
* @param {String} fieldId the ID of the field this widget is being used for
**/
export function mmpKetcherHandler(fieldId) {
  const frameGroup = document.getElementById(`${fieldId}-frame-group`)
  const loadingMsg = document.getElementById(`${fieldId}-loading`);
  const ketcherFrame = document.getElementById(`${fieldId}-frame`);
  if (ketcherFrame.contentWindow.ketcher) {
    // Hide the loading message and show the iframe and highlight buttons
    loadingMsg.style.display = "none"
    frameGroup.style.display = "";
    ketcherLoaded(fieldId);
  } else {
    // iframe not loaded, check status again after 50 milliseconds
    window.setTimeout(function () { mmpKetcherHandler(fieldId) }, 50);
  }
}
/**
 * Called when the ketcher iframe has loaded, this function:
 * - retrieves the ketcher object
 * - initializes widgetData and populates ketcher with data from the previous form post when applicable
 * - starts event listeners to keep widgetData up to date and to manage variable region highlights
 * @param {String} fieldId the ID of the field this widget is being used for
 */
async function ketcherLoaded(fieldId) {
  const ketcherFrame = document.getElementById(`${fieldId}-frame`)
  const ketcher = ketcherFrame.contentWindow.ketcher;
  const textarea = document.getElementById(fieldId);
  var widgetData = {
    moltext: "",
    variable: {atoms: [], bonds: []},
    maps: { atoms: {}, bonds: {}}
  };
  if (textarea.value) {
    widgetData = JSON.parse(textarea.value);
    widgetData = refreshWidgetData(widgetData)
    await ketcher.setMolecule(widgetData.moltext);
    // Highlight variable region (if highlighted in widgetData)
    try{
      ketcher.editor.selection({ atoms: widgetData.variable.atoms, bonds: widgetData.variable.bonds });
      updateHighlights('variable', ketcher, widgetData, textarea)
    }catch(error){
      // If Ketcher throws an error when trying to highlight atoms and bonds (as sometimes occurs with
      // enantiomers), clear the saved variable region.
      widgetData.variable.atoms = [];
      widgetData.variable.bonds = [];
    }
    
  }
  ketcher.editor.subscribe("change", function () {
    onKetcherChange(ketcher, widgetData, textarea)
  });
  document.getElementById(`${fieldId}-variable-btn`).addEventListener("click", function () {
    updateHighlights('variable', ketcher, widgetData, textarea);
  });
  document.getElementById(`${fieldId}-clear-btn`).addEventListener("click", function () {
    clearHighlights(ketcher, widgetData, textarea);
  });
}

/**
 * Called every time the molecule in the ketcher widget changes, this function updates widgetData.moltext
 * with the newly drawn structure & removes deleted atoms & bonds from widgetData.variable
 * @param {Object} ketcher the ketcher object
 * @param {Object} widgetData the widgetData that contains the atoms and bonds used in the variable region
 * @param {Object} textarea the textarea for this widget that contains widgetData
 */
function onKetcherChange(ketcher, widgetData, textarea){
  ketcher.getMolfile().then(moltext => {
    widgetData.moltext = moltext
    const currentMol = ketcher.editor.selection("all")
    ketcher.editor.selection(null)
    widgetData.variable.atoms = widgetData.variable.atoms.filter(atom => currentMol.atoms.includes(atom))
    widgetData.variable.bonds = widgetData.variable.bonds.filter(bond => currentMol.bonds.includes(bond))
    textarea.value = JSON.stringify(widgetData)
  })
}

/**
 * Called when the variable button is clicked, this retrieves the currently
 * highlighted substructure from the ketcher widget, highlights it in the corresponding color,
 * and updates `widgetData`. 
 * @param {String} region the region being highlighted ("constant" or "variable")
 * @param {Object} ketcher the ketcher object
 * @param {Object} widgetData the widgetData that contains the atoms and bonds used in the variable & constant regions
 * @param {Object} textarea the textarea for this widget that contains widgetData
 */
function updateHighlights(region, ketcher, widgetData, textarea) {
  const ketcherSelection = ketcher.editor.selection();
  if (!ketcherSelection) {
    // If nothing is selected, don't update
    return
  }
  updateMaps(ketcher, widgetData, textarea);
  // Add the currently selected atoms and bonds to the widgetData
  ["atoms", "bonds"].forEach(componentType => {
    const components = ketcherSelection[componentType] || []
    components.forEach(component => {
      if (!widgetData[region][componentType].includes(component)) {
        widgetData[region][componentType].push(component)
      }
    })
  })
  // Update the value of the textarea & highlight the molecule
  textarea.value = JSON.stringify(widgetData)
  const color = (region == "constant") ? '#E2E2E2' : '#f59043'
  const highlight = { atoms: widgetData[region].atoms, bonds: widgetData[region].bonds, color: color }
  ketcher.editor.highlights.create(highlight);
}

/**
 * Clears any highlights made by the user and updates the widgetData accordingly
 * @param {Object} ketcher the ketcher object
 * @param {Object} widgetData the widgetData that contains the atoms and bonds used in the variable region
 * @param {Object} textarea the textarea for this widget that contains widgetData
 */
function clearHighlights(ketcher, widgetData, textarea) {
  ketcher.editor.highlights.clear();
  widgetData.variable.atoms = [];
  widgetData.variable.bonds = [];
  textarea.value = JSON.stringify(widgetData)
}


/**
 * Ketcher atom/bond indices do not directly correspond to the index of the related atom/bond in moltext
 * because when components of molecules are deleted or rotated, ketcher continues incrementing the indices
 * from the last used number (Example: draw benzene, atoms have indices 0-5, delete benzene and redraw, atoms have indices 6-11).
 * As a result, we need to keep track of which ketcher atom/bond indices map to which atom/bond indices in the final moltext
 * (this is accomplished in widgetData.maps.atoms and widgetData.maps.bonds).
 * 
 * This function updates widgetData.maps with the current mappings of {ketcherIndex: moltextIndex}.
 * @param {Object} ketcher the ketcher object
 * @param {Object} widgetData the widgetData that contains the atom and bond index maps
 * @param {Object} textarea the textarea for this widget that contains widgetData
 */
function updateMaps(ketcher, widgetData, textarea) {
  const ketcherSelection = ketcher.editor.selection('all');
  ketcher.editor.selection(null);
  ["atoms", "bonds"].forEach(componentType => {
    if (ketcherSelection[componentType] !== undefined) {
      for (var i = 0; i < ketcherSelection[componentType].length; i++) {
        widgetData.maps[componentType][ketcherSelection[componentType][i]] = i ;
      }
    }
  })
  textarea.value = JSON.stringify(widgetData)
}

/**
 * This function is called when ketcher loads if there is already widgetData (from the previous form post).
 * In the previous form post, the atom & bond numbers in the sketcher likely don't start from 0, incrementing by 1 with
 * each additional atom/bond, because when the user deletes & rotates pieces in the sketcher, the atom numbers change. This
 * is why we keep atom & bond index maps in widgetData.
 * 
 * As a result, the constant and variable data in widgetData from the previous form post might include atom/bond indices that are greater
 * than the total number of atoms/bonds in the molecule. Because the page reloaded with form errors, the ketcher instance has
 * been refreshed, starting its atom/bond indices from 0, incrementing by 1 with each added atom/bond, so these greater indices
 * do not exist anymore. To address this, we use the existing atom/bond index maps to update the constant & variable data in widgetData
 * with new indices that correspond to the correct atoms in the new ketcher instance.
 * 
 * Example:
 * - A user draws benzene, deletes it, redraws benzene, highlights the entire molecule as the constant region
 * and submits the form. The molecule has 6 atoms and 6 bonds each numbered 6-11 (because 0-5 were used by the
 * original benzene that was deleted). widgetData.constant.atoms = [6,7,8,9,10,11] & widgetData.maps.atoms = {6:0, 7:1, 8:2, 9:3, 10:4, 11:5}
 * - The form returns with errors because no variable region was included. The molecule is redrawn in ketcher,
 * but now its atom/bond indices are 1-5 (because the atom/bond index counts start over from 0). widgetData.constant.atoms
 * still has the previous value [6,7,8,9,10,11], but it should be [0,1,2,3,4,5].
 * @param {Object} widgetData a dictionary of widgetData 
 */
function refreshWidgetData(widgetData){
  ["atoms", "bonds"].forEach(componentType => {
    widgetData["variable"][componentType] = widgetData["variable"][componentType].map(function (component){
      return widgetData.maps[componentType][component]
    })
  })
  return widgetData
}



