
/**
 * This function should be called in the 'row-reordered' event trigger on a dataTable with
 * rowReorder = true. It saves the new user-specified order to the database so that it is 
 * persistent across pages/reloads. For the propcalc table, this function also re-orders
 * the compounds in the grid to match the order in the table.
 * @param {String} tableId the id of the table element whose CompoundOccurrence rows can be reordered
 */
function saveCoReordering(tableId){
    const table = document.getElementById(tableId)
    const coOrder = []
    for (var i = 1; i < table.rows.length; i++) {
        coOrder.push(parseInt(table.rows[i].id.slice(3)))
    }
    // Save order to DB
    const collectionPk = JSON.parse(document.getElementById('collection-pk').textContent);
    $.ajax({
        url: `/ajax/update-collection/${collectionPk}/`,
        type : "post",
        data: JSON.stringify({"update_field": "co_order", "co_order": coOrder}),
    })
    if(tableId == "propcalc-table"){
        // Update order of grid to match table order
        const gridContainer = document.getElementById("propcalc-grid-container")
        for(var i = coOrder.length-1; i > 0; i--){
            const beforeElement = document.querySelector(`.propcalc-grid-item#co-${coOrder[i-1]}`)
            const afterElement = document.querySelector(`.propcalc-grid-item#co-${coOrder[i]}`) 
            gridContainer.insertBefore(beforeElement, afterElement)
        }
    }
}

/**
 * This function should be called in the `order.dt` event trigger on a dataTable with
 * rowReorder = true. It disables row reordering when the the table is not sorted by the order
 * column (and enables it when it is sorted exclusively by the order column).
 * @param {String} tableId the id of the table element whose CompoundOccurrence rows can be reordered
 */
function toggleRowReorder(tableId){
    const dataTable = $(`#${tableId}`).DataTable()
    var order = dataTable.order();
    // If the table is ordered exclusively by the "order" column (assumed to be column 0)
    if((order[0] == 0) ||(order.length == 1 && order[0][0] == 0)){
        // Allow users to reorder the rows using drag & drop
        dataTable.rowReorder.enable();
        // Remove the 'reorder' class from the first column so the cursor doesn't suggest that the row can be moved
        $(`table#${tableId} td.order-column`).toggleClass("reorder", true)
    }else{
        dataTable.rowReorder.disable();
        $(`table#${tableId} td.order-column`).toggleClass("reorder", false)
    }
}

/**
 * This function should be called after DataTable initialization for all DataTables with rowReorder: true.
 * It adds two event listeners to the table. The first saves the user specified order to the database every
 * time a row is moved. The second disables reordering when the table is not sorted by the
 * order column, since that leads to confusing behavior in the UI
 * @param {String} tableId the id of the table element whose CompoundOccurrence rows can be reordered
 */
export function enableRowReorder(tableId){
    const dataTable = $(`#${tableId}`).DataTable()
    dataTable.on('row-reordered', function () {saveCoReordering(tableId)});
    dataTable.on('order.dt', function () {toggleRowReorder(tableId)});
}