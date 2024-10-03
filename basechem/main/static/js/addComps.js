import { clearFormMessages } from "../../../mni_common/static/mni_common/js/ajaxForms.js"
/**
 * This is called when the AddCompoundForm is submitted and the ajax form post returns. It adds
 * rows to the analysis table with the newly added compounds and resets the ketcher widget that
 * was used to add the new compounds.
 * @param {String} analysis the name of the current analysis (ex: "align")
 * @param {Array} tableRows an array of strings containing HTML for `tr` elements
 * @param {Array} gridItems an array of strings containing HTML for grid elements (propcalc only)
 */
export function addComps(analysis, tableRows, gridItems){
  // Add new rows to the DOM
  const table = document.getElementById(`${analysis}-table`)
  const dataTable = $(`#${table.id}`).DataTable()
  var orderNum = parseInt(dataTable.column(0).data().sort().reverse()[0]) + 1
  tableRows.forEach(tableRow => {
    const row = dataTable.row.add($(tableRow))
    dataTable.cell(row.index(), 0).data(orderNum)
    orderNum += 1
  })
  dataTable.draw()
  // Add new grid items to the grid (propcalc only)
  gridItems.forEach(gridItem => {
    const grid = document.getElementById('propcalc-grid-container')
    grid.insertAdjacentHTML("beforeEnd", gridItem)
  })
  // Clear the sketcher
  const form = document.getElementById("add-comp-form")
  const ketcher = form.querySelector("iframe").contentWindow.ketcher
  ketcher.setMolecule("")
}

/**
 * Adds event listeners for adding compounds to collections
 */
export function addCompsEventListeners(){
  $(".add-comp-btn").on("click", addCompOnClick)
  $('#add-comp-modal').on('hidden.bs.modal', function () {
    clearFormMessages(document.getElementById("add-comp-form"))
  })

}

/**
 * Handles opening the add-comp-modal when an add-comp-button is clicked. This is used instead
 * of a basic bootstrap HTML modal toggle so that the ketcher widget in the modal can be pre-populated
 * with moltext when appropriate.
 */
export function addCompOnClick(){
  const form = document.getElementById("add-comp-form")
  const ketcher = form.querySelector("iframe").contentWindow.ketcher
  const moltext = this.getAttribute("data-moltext")
  try{
    ketcher.setMolecule(moltext)
  }catch(error){
    if(moltext){
      alert("Please wait a few seconds and try again - the sketcher is still loading.")
      return
    }
  }
  $("#add-comp-modal").modal("show")
}

/**
 * This listens for new rows added to the results table (specified w/ tableId). This happens
 * when new molecules are added to the collection. When a new row is added, calls `helperFunc` to handle
 * any analysis-specific side effects, including collecting task results for the new compound.
 * @param {String} tableId the id of the table element where new compounds are added
 * @param {Object} parsedData an object with all results for the current analysis
 * @param {Object} viewer the displayed 3dmol viewer
 * @param {Object} helperFunc an analysis-specific helper function that should accept trElement, parsedData, and viewer
 *  as arguments. This function should handle anything that happens when a new compound is added to a collection.
 */
export function handleNewCompsAdded(tableId, parsedData, viewer, helperFunc){
  const tBody = document.getElementById(tableId).querySelector("tbody");
  // Callback function to execute when mutations are observed
  const callback = (mutationList, observer) => {
    for (const mutation of mutationList) {
      mutation.addedNodes.forEach(node => {
        if((node.tagName=="TR") && (node.getAttribute("data-new-row") == "true" )){
          node.setAttribute("data-new-row", "false")
          helperFunc(node, parsedData, viewer)
          $(node.querySelector(".add-comp-btn")).on("click", addCompOnClick)
          $(".hover-tooltip").tooltip({trigger:"hover"})
        }
      })
    }
  };
  const observer = new MutationObserver(callback);
  observer.observe(tBody, {childList: true});
}