const PRIMARY_OPACITY = 1;
const SECONDARY_OPACITY = 0.7;
const STYLE = 'stick';
const CARBON_COLORS = [0xFFE8BC, 0xB3FFE0, 0xFFCEEA, 0x87F07D, 0x8CA7AB, 0x88D9FF, 0xB1AE70, 0x8882AD]
const ATOM_COLORS = {
    'O': 0xFE2103,
    'S': 0xFFE806,
    'Cl': 0X32A852,
    'F': 0x1EF01F,
    'N': 0x2333FF,
}
// All torsions alerts <= TORSION_ALERTS_GREEN_CUTOFF will be colored green
const TORSION_ALERTS_GREEN_CUTOFF = 0.5
// All torsions alerts >= TORSION_ALERTS_GREEN_CUTOFF and <= TORSION_ALERTS_ORANGE_CUTOFF
// will be colored orange. Anything greater than this value will be colored red.
const TORSION_ALERTS_ORANGE_CUTOFF = 1.5


/***************************************************
 *                 Manage Models                   *
 ***************************************************/

/**
 * Removes a model from the viewer
 * @param {Object} models a dictionary mapping model ids to their 3dmol models
 * @param {Object} colorMap a dictionary mapping model ids to the color used for that model
 * @param {Object} colorFrequency a dictionary mapping color hex codes to the number
 *  of times that color has been used in the current viewer
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @param {String} id the basechem identifier of the model being removed
 * @param {Boolean} deleteUnusedColor if true, removes the color of the deleted model from the colorMap
 *      (this should usually be true, if it is false, make sure to update the colorMap somewhere else)
 */
 function removeModel(models, colorMap, colorFrequency, viewer, id, deleteUnusedColor=true){
    // Remove existing model, torsionShapes, and toklatShapes from viewer and delete from models dict
    models[id].torsionShapes.forEach(shape => { viewer.removeShape(shape)})
    models[id].toklatShapes.forEach(shape => { viewer.removeShape(shape)})
    viewer.removeModel(models[id].model)
    delete models[id];
    // Delete color from color map and re-color other tags
    if(deleteUnusedColor){
        colorFrequency[colorMap[id]] -= 1
        delete colorMap[id]
    }
    colorTags(colorMap)
}


/**
 * Adds a styled model to the viewer
 * @param {Object} models a dictionary mapping model ids to their 3dmol models
 * @param {Object} colorMap a dictionary mapping model ids to the color used for that model
 * @param {Object} colorFrequency a dictionary mapping color hex codes to the number
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @param {String} moltext the moltext of the model to add
 * @param {String} id the basechem identifier of this model (ex s-3 for series or 121-0 for a conformer)
 * @param {String} format "sdf" or "pdb", what is the format of the moltext
 * @param {Float} opacity the opacity of the model
 * @param {String} style "cartoon" or "stick", the atom style of the model
 * @param {Object} torsionAlerts a dictionary of torsion alert information that maps a pair
 *  of atom indices to the energy (ex: {"6,7": "5.9"})
 * @returns the model added to the viewer
 */
function addModel(models, colorMap, colorFrequency, viewer, moltext, id, {format="sdf", opacity=PRIMARY_OPACITY, style=STYLE, torsionAlerts ={}, toklatAnnotations={}} = {}){
    // Create new model
    const m = viewer.addModel(moltext,format);
    // Remove existing model, torsionShapes, and toklatShapes if they exist
    if(models[id]){
        viewer.removeModel(models[id].model)
        models[id].torsionShapes.forEach(shape => { viewer.removeShape(shape)})
        models[id].toklatShapes.forEach(shape => { viewer.removeShape(shape)})
    }
    // Add model to dictionary and style model
    models[id] = {"model": m, "torsionShapes": [], "toklatShapes": []}
    styleModel(models, colorMap, colorFrequency,id, {style:style, opacity:opacity})
    addTorsionAlerts(viewer, models, id, torsionAlerts)
    addToklatAnnotations(viewer, models, id, toklatAnnotations)
    return m;
}

/***************************************************
 *                   Style Models                  *
 ***************************************************/

/**
 * Adds the styling elements to a 3dmol model
 * @param {Object} models a dictionary mapping model ids to their 3dmol models
 * @param {Object} colorMap a dictionary mapping model ids to the color used for that model
 * @param {Object} colorFrequency a dictionary mapping color hex codes to the number
 * @param {String} id the basechem identifier of the model being styled
 * @param {String} style "cartoon" or "stick"
 * @param {Float} opacity_val the opacity of the model
 * @param {Object} atomColorDict a dictionary mapping periodic table elements to the hex code that should be used
 */
 function styleModel(models, colorMap, colorFrequency, id, {style=STYLE, opacity_val=PRIMARY_OPACITY, atomColorDict=ATOM_COLORS} = {}){
    const m = models[id].model
    var options = {}
    if(style=="cartoon"){
        options ={cartoon:{opacity: opacity_val}};
    }else if(style=="stick"){
        options = {stick:{opacity: opacity_val}};
    }
    // Apply style
    m.setStyle({}, options);
    // Hide all hydrogens
    const hs = {elem: ["H", "Hsh","Hlp","Hh","Hp", "Hl", "Hs"]}
    m.setStyle(hs, {stick:{hidden: true}});
    // Set atom colors
    const atomColors = JSON.parse(JSON.stringify(atomColorDict));
    if(id != "receptor"){
        atomColors['C'] = pickCarbonColor(colorMap, colorFrequency, id);
        colorTags(colorMap) // Update the select2 tags to have the appropriate background colors
    }
    m.setColorByElement({}, atomColors);
}

/**
 * Determine if torsion alerts are currently displayed in the viewer
 * @param {Object} models a dictionary mapping model ids to their 3dmol models
 * @returns a boolean, are torsion alerts currently shown in the viewer?
 */
function areTorsionAlertsShown(models){
    var showAlerts = false
    for(const mId in models){
        models[mId].torsionShapes.forEach(shape => {
            if(!shape.hidden){
                showAlerts = true
            }
        })
    }
    return showAlerts
}

/**
 * Determine if toklat annotations are currently displayed in the viewer
 * @param {Object} models a dictionary mapping model ids to their 3dmol models
 * @returns a boolean, are toklat annotations currently shown in the viewer?
 */
function areToklatAnnotationsShown(models){
    var showAnnotations = false
    for(const mId in models){
        models[mId].toklatShapes.forEach(shape => {
            if(!shape.hidden){
                showAnnotations = true
            }
        })
    }
    return showAnnotations
}

/**
 * Add torsion alerts (in the form of colored cylinders) to a model in the viewer.
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @param {Object} models a dictionary mapping model ids to their 3dmol models
 * @param {String} id the basechem identifier of the model being styled
 * @param {Object} alerts a dictionary of torsion alert information that maps a pair of
 *  atom indices (ex: "6,7") to the energy (ex: "5.9")
 */
function addTorsionAlerts(viewer, models, id, alerts){
    const m = models[id].model
    for(var [atomIndices, energy] of Object.entries(alerts)){
        atomIndices = atomIndices.split(",")
        energy = parseFloat(energy)
        const showAlerts = areTorsionAlertsShown(models)
        const atom1 = m.selectedAtoms({serial:atomIndices[0]})[0]
        const atom2 = m.selectedAtoms({serial:atomIndices[1]})[0]
        var color = "red"
        if(energy <= TORSION_ALERTS_ORANGE_CUTOFF){color = "orange"}
        if(energy <= TORSION_ALERTS_GREEN_CUTOFF){color = "green"}
        const cylinder = viewer.addCylinder({
            start:{x:atom1.x, y:atom1.y , z:atom1.z},
            end:{x:atom2.x ,y:atom2.y ,z:atom2.z},
            radius:0.4,
            fromCap:1,
            toCap:1,
            color:color,
            opacity: 0.7,
            hidden: !showAlerts,
        })
        models[id].torsionShapes.push(cylinder)
    }
}

/**
 * Add Toklat annotations (in the form of colored cylinders and wireframe spheres) to a model in the viewer.
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @param {Object} models a dictionary mapping model ids to their 3dmol models
 * @param {String} id the basechem identifier of the model being styled
 * @param {Object} annotation a dictionary of Toklat annotation data of the form {"cylinders": [{"start", "end", "color", "radius"}], "spheres": [{"center"}]}
 */
function addToklatAnnotations(viewer, models, id, annotations){
    if(Object.keys(annotations).length == 0){
        return
    }
    const showAnnotations = areToklatAnnotationsShown(models)
    annotations.cylinders.forEach(cylinderDict => {
        const cylinder = viewer.addCylinder({
            start:cylinderDict.start,
            end:cylinderDict.end,
            color:cylinderDict.color,
            radius:cylinderDict.radius,
            fromCap:1,
            toCap:1,
            opacity: 0.7,
            hidden: !showAnnotations,
            dashed: true,
        })
        models[id].toklatShapes.push(cylinder)
    })

    annotations.spheres.forEach(sphereDict => {
        const sphere = viewer.addSphere({
            center: sphereDict.center,
            color: "grey",
            radius: 0.5,
            hidden: !showAnnotations,
            wireframe: true,
        })
        models[id].toklatShapes.push(sphere)
    })

}

/**
 * Color the select2 tags so they match the color of the carbons for each selection. The tag
 * DOM elements are reset every time a new selection is made, so this function re-applies
 * the color to every tag each time it is called (instead of just adding color to the new tag)
 * @param {Object} colorMap a dictionary mapping model ids to the color used for that model
 */
function colorTags(colorMap){
    for(var id in colorMap){
        const hex_color = parseInt(colorMap[id]).toString(16)
        var tag = document.querySelector(`.select2-selection__choice[data-id='${id}']`);
        if(!tag){
            tag = document.querySelector(`button.toggle-compound#${id}`)
        }
        if(tag){
            tag.style.backgroundColor = "#"+hex_color;
        }
    }
}

/**
 * Sorting function that sorts an array of [[color_key, color_frequency]] by frequency,
 * with the least frequent colors first. Ties are broken by the order that the color
 * appears in the `CARBON_COLORS`
 * @param {Array} a [color_key, color_frequency] of the first color
 * @param {Array} b [color_key, color_frequency] of the second color
 * @returns an integer, the sorting value to use
 */
function sortColors(a,b){
    var sortVal = a[1]-b[1]
	if(a[1]-b[1] == 0){
  	    sortVal = CARBON_COLORS.indexOf(parseInt(a[0])) - CARBON_COLORS.indexOf(parseInt(b[0]))
    }
    return sortVal
}


/**
 * Styles the receptor model so that there is more detail in the binding pocket
 * @param {Object} models a dictionary mapping model ids to their 3dmol models
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @param {Array} residues an array of binding pocket residues (integers)
 * @param {Number} slabNear minimum value for the viewer's slab
 * @param {Number} slabFar maximum value for the viewer's slab
 */
function styleReceptor(models, viewer, residues, slabNear, slabFar){
    const receptor = models["receptor"].model
    const pocket = {rescode:residues}
    receptor.setStyle({}, {cartoon:{opacity: 0.5}})
    receptor.setStyle(pocket, {stick:{opacity: 0.9}, cartoon:{opacity: 0.5}});
    viewer.setSlab(slabNear,slabFar)
    viewer.render()
}

/**
 * Chooses the least used color from `CARBON_COLORS` for `id`'s model. Also updates the
 * color that is being used for `id`'s model in the `colorMap`
 * @param {Object} colorMap a dictionary mapping model ids to the color used for that model
 * @param {Object} colorFrequency a dictionary mapping color hex codes to the number
 * @param {String} id the basechem identifier of the model being colored
 * @param {Boolean} hex if true, returns the color as a hex string (#ffffff). If false,
 *      returns as an integer
 * @returns the color as either an integer or hex string
 */
function pickCarbonColor(colorMap, colorFrequency, id, hex=false){
    const colorId = getColorId(id)
    let color = null
    if(colorId in colorMap){
        color = colorMap[colorId]
    }else{
        const sortedColors = Object.entries(colorFrequency).sort((a, b) => sortColors(a,b));
        color = sortedColors[0][0]
        colorFrequency[color] += 1
        colorMap[colorId] = color 
    }
    if(hex){
        return `#${parseInt(color).toString(16)}`
    }else{
        return color
    }
}

/**
 * This function returns the id (string) that should be used to look up the color in the colorMap.
 * Most of the time the returned id is the same as the original, but in the case of torsion scans
 * the color is determined by the compound occurrence, not each model individually.
 * @param {String} id  the basechem identifier of a model/conformer (ex: "co-7" or "c2-co7-12")
 * @returns a string, the id that should be used to look up the color in the colorMap
 */
 function getColorId(id){
    let colorId = id
    const groupName = JSON.parse(document.getElementById('group-name').textContent)
    if(groupName.includes("torsion")){
        // Torsion uses the same color for all models from the same compound occurrence
        var matches = id.match(/co-\d+/)
        if(matches){
            colorId = matches[0]
        }else{
            matches = id.match(/co\d+/)
            colorId = matches[0].slice(0, 2) + "-" + matches[0].slice(2);
        }  
    }
    return colorId
}



/***************************************************
 *                    Downloads                    *
 ***************************************************/

/**
 * Adds event listeners to the download buttons
 */

/**
 * Adds event listeners to download buttons. Requires that elements with `querySelector`
 * have a `data-url` attribute with the download url endpoint. If the `id` of the element
 * is 'displayed', appends (a_b..._z) to the end of the url, where a,b,...z are the 
 * basechem ids of the models that are currently being displayed (excluding receptors)
 * @param {Object} models a dictionary mapping model ids to their 3dmol models
 * @param {String} querySelector a querySelector that is unique to all the
 *  download buttons (ex. '.dock-download')
 */
function download(models, querySelector){
    $(querySelector).on('click', function(){
        var url = this.getAttribute('data-url')
        if(this.id == "displayed"){
            // Add the ids of the displayed models to the end of the url
            var displayed = []
            for(var key in models){
                if(key != "receptor"){
                    displayed.push(key)
                }
            }
            url = url + displayed.join("_")
        }
        else if (this.id == "chart"){   
            var chartIds = []
            const chart = Chart.getChart('dihedral-energy-chart');
            for(let datasetIndex=1; datasetIndex<chart.data.datasets.length; datasetIndex++){
                const dataset = chart.data.datasets[datasetIndex]
                var coId = dataset.id; 
                for(let pointIndex=1; pointIndex<dataset.data.length; pointIndex++){
                    var dihedral = dataset.data[pointIndex].x
                    chartIds.push(getTorsionConfId(dataset, dihedral))
                }
                
            }
            url = url + chartIds.join("_")
        }
        window.location = url;
    })
}

/***************************************************
 *                      Util                       *
 ***************************************************/

/**
 * Zooms the given viewer to the coordinate space of the last model in `models` when a button
 * (specified by `querySelector`) is clicked.
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @param {Object} models a dictionary mapping model ids to their 3dmol models
 * @param {Object} querySelector a querySelector that is unique to the fit-viewer button
 * @param {Number} slabNear (optional) near boundary of the slab (ex: -10)
 * @param {Number} slabFar (optional) far boundary of the slab (ex: 10)
 */
function fitViewerToModels(viewer, models, querySelector, slabNear, slabFar){
    $(querySelector).on('click', function(){
        // Pick the last model in `models`
        const mIds = Object.keys(models)
        let mId = mIds.pop()
        if(mId == "receptor"){
            // Don't fit to the receptor, pick the next model
            mId = mIds.pop()
        }
        const m = models[mId].model
        viewer.zoomTo({model: m})
        if(slabNear & slabFar){
            viewer.setSlab(slabNear, slabFar)
        }
        viewer.render()
    })
}

/**
 * Adds an event listener to the button matching the query selector that will toggle torsion alerts on/off.
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @param {Object} models a dictionary mapping model ids to their 3dmol models
 * @param {String} querySelector the query selector for the button toggling torsion alerts
 */
function toggleTorsionAlerts(viewer, models, querySelector){
    $(querySelector).on('click', function(){
        const showAlerts = !areTorsionAlertsShown(models)
        for(const mId in models){
            models[mId].torsionShapes.forEach(shape => {
                shape.updateStyle({hidden: !showAlerts})
            })
        }
        viewer.render()
    })
}

/**
 * Adds an event listener to the button matching the query selector that will toggle Toklat annotations on/off.
 * @param {Object} viewer the 3dmol viewer that is displayed
 * @param {Object} models a dictionary mapping model ids to their 3dmol models
 * @param {String} querySelector the query selector for the button toggling Toklat annotations
 */
function toggleToklatAnnotations(viewer, models, querySelector){
    $(querySelector).on('click', function(){
        const showAlerts = !areToklatAnnotationsShown(models)
        for(const mId in models){
            models[mId].toklatShapes.forEach(shape => {
                shape.updateStyle({hidden: !showAlerts})
            })
        }
        viewer.render()
    })
}

/**
 * Given a dataset and dihedral, return the expected conformer ID
 * @param {Object} dataset a dataset from a ChartJS chart
 * @param {Number} dihedral an number between -180 and 180 (inclusive)
 * @returns a string, the expected ID of the 3dmoljs model in the viewer
 */
function getTorsionConfId(dataset, dihedral){
    const prefix = document.querySelector(`tr.compound-row#${dataset.id}`).getAttribute("data-conf-prefix")
    var modelId = `${prefix}-${dihedral}`
    if(dihedral < 0){
        modelId = `${prefix}-neg${Math.abs(dihedral)}`
    }
    return modelId
}


export {PRIMARY_OPACITY, SECONDARY_OPACITY, STYLE, CARBON_COLORS, ATOM_COLORS, download, removeModel, addModel, colorTags, fitViewerToModels, toggleTorsionAlerts, toggleToklatAnnotations, pickCarbonColor, styleReceptor, getTorsionConfId}
