/**
   * This is called after the save gems form is submitted and the ajax form post returns. It
   * updates the options in the conformer dropdown in the align/dock table to reflect the currently saved
   * gems
   * @param {Integer} coPk PK of compound occurrence whose conformers need to be updated
   * @param {Array} savedConfs array of conformer IDs for the currently saved gems (ex: ["co-512-1", "co-512-29"])
   * @param {Array} unsavedGems array of conformer IDs for previously save gems that have been un-saved (ex: ["co-512-5"])
   * @param {String} selectElementClass the class name of the select HTML element whose options need to be updated
   */
 export function updateSavedGems(cPk, coPk, savedConfs, unsavedGems, selectElementClass){
    hideUnsavedGems(unsavedGems)
    // Make sure the checkboxes in the save gems forms are checked
    savedConfs.forEach(confId => {
      var checkbox = document.querySelector(`input[type=checkbox][name=gems][value="${confId}"]`)
      if(!checkbox){
        checkbox = document.querySelector(`input[type=checkbox][name=other_gems][value="${confId}"]`)
      }
      checkbox.checked = true
    })
    if(!selectElementClass){
      return
    }
    const options = document.querySelectorAll(`option[data-id^="c${cPk}-co${coPk}"]`)
    const selectElement = document.querySelector(`.${selectElementClass}#co-${coPk}`)
    for(let i=0; i<options.length; i++){
        const op = options[i]
        var updatedInnerText = op.innerText
        if(savedConfs.includes(op.getAttribute("data-id")) && !op.innerText.includes("*")){
            // Conformer is newly saved, add asterisks to dropdown
            updatedInnerText = op.innerText +" *"
        }else if(!savedConfs.includes(op.getAttribute("data-id")) && op.innerText.includes("*")){
            // Conformer was previously saved, but is not saved anymore: remove asterisks from dropdown
            updatedInnerText = op.innerText.substring(0, op.innerText.length - 2)
        }
        op.innerText = updatedInnerText
        // Make sure changes are reflected in the select2 dropdown
        const dataId = op.getAttribute('data-id')
        $(selectElement).find(`[data-id="${dataId}"]`).data('display', updatedInnerText);
    }
    $(selectElement).trigger('change');
  }


  /**
   * This is called from updateSavedGems to remove options from the "Previously saved" form field
   * when previously saved gems are unsaved
   * @param {Array} gemIds array of gem IDs for previously saved gems that have been un-saved (ex: ["co-512-5"])
   */
  function hideUnsavedGems(gemIds){
    gemIds.forEach(gemId => {
        const checkbox = document.querySelector(`input[name='other_gems'][value='${gemId}']`)
        if(checkbox){
          checkbox.parentElement.style.display = "none"
        }
    })
  }

  export function saveViewerGems(models, querySelector){
    $(querySelector).on("click", function(){
      // Get the groupName of the tasks that contain the moltext for the selected conformers
      const groupName = JSON.parse(document.getElementById('group-name').textContent)
      // Get array of all conformers in the viewer, excluding receptors and series
      var viewerConfIds = Object.keys(models).filter(function(confId) {
        const regex = new RegExp('^(s-\\d+)|(receptor)$');
        return !regex.test(confId);
      });
      $.ajax({
        url : $(this)[0].getAttribute("data-post-url"),
        type : "post",
        data : JSON.stringify({
          "groupName": groupName,
          "viewerConfIds": viewerConfIds
        }),
        processData: false,
        contentType: false,
        success : function(response) {
          if(!response.errors.length){
            // Success
            response.success_methods.forEach(method => eval(method))
          }else{
            // Known failure
            alert(response.errors.join("\n"))
          }
        },error: function(){
          // Unknown failure
          alert("An unexpected error occurred while saving gems, contact admins.")
        }
      });
    })
  }