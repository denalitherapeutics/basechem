/**
 * Adds a row below the headings to be used as filters
 * @param {string} tableId the id of the table element
 */
export function addFilterRow(tableId){
    const filterRow = $(`#${tableId} thead tr`)
        .clone(true)
        .addClass('filters')
        .appendTo(`#${tableId} thead`);
    filterRow[0].id = `${tableId}-filters`
}

/**
 * This adds filtering inputs to a DataTable. This should be called within the initComplete
 * function. The filtering inputs are styled by and rely on select2.
 * @param {Array} blankFilter array of strings where each string is the name of a table
 *      column whose filter search bar should not pre-populate with options.
 *      One reason this may be used for columns where each value is expected to be unique
 * @param {Array} hideFilter an array of strings where each string is the name of a
 *      table column that should not have a filter input field
 * @param {Array} exactMatch an array of strings where each string is the name of a
 *      table column whose query should be an exact match for its results (ex: number columns)
 */
export function filterTable({ blankFilter =[], hideFilter=[], exactMatch=[] } = {}){
    const tableId = this[0].id
    const api = this.api();
    // For each column in the datatable
    api.columns().eq(0).each(function (colIdx) {
        const filterCell = $(`#${tableId}-filters th`).eq($(api.column(colIdx).header()).index());
        const title = $(filterCell).text();
        // Use the column title to create an HTML ID for the select field.
        // We replace spaces and unusual characters with the empty string
        var elementId = title.toLowerCase().replace(/\s|\(|\)|\'|,|\//g, "").replace('#', "num");
        if(hideFilter.includes(elementId)){
            // This column should not have an input field, stop here
            $(filterCell)[0].innerHTML = ""
            return;
        }
        // Generate the dropdown options based on the unique values in the column
        var options = ""
        if(!blankFilter.includes(elementId)){
            api.column(colIdx).data().unique().sort().each( function ( d, j ) {
                if(!d.includes("\n")){
                    options = options + `<option>${d}</option>`
                }
            })
        }
        // Add a select widget to the filterCell
        $(filterCell).html(`<select id="${tableId}_${elementId}_select"  data-exact-match=${exactMatch.includes(elementId)} class="table-filter" multiple>${options}</select>`);
        // Style using select2
        $(`select#${tableId}_${elementId}_select`).select2({
            placeholder:title,
            allowClear:true,
            tags:true,
            width:"100%",
            // Override "No results" with the empty string
            language: {"noResults": function (){return ""}},
            // Sort the options in alphabetical order
            sorter: data => data.sort((a, b) => a.text.localeCompare(b.text, navigator.languages[0] || navigator.language, {numeric: true, ignorePunctuation: true})),
        });
        // Hide the dropdown if there are no options to choose from
        if(!options){
            $(`select#${tableId}_${elementId}_select`).on('select2:opening select2:close', function(e){
            $('body').toggleClass('select2-no-options', e.type=='select2:opening');
            });
        }
        // On change, filter the table
        $('select',$(`#${tableId}-filters th`).eq($(api.column(colIdx).header()).index()))
            .on('change', function (e) {
                e.stopPropagation();
                const selectElement = document.getElementById(`${tableId}_${elementId}_select`)
                // Search the column for values that match the selected values
                const selectedValues = Array.from(this.selectedOptions).map(({ value }) => value.replace(/[-[\]{}()*+?.,\\^$|#\s]/g, '\\$&'));
                var regex = this.value ? "("+selectedValues.join("|")+")" : ''
                if(selectElement.getAttribute("data-exact-match") === "true"){
                    regex = this.value ? "^("+selectedValues.join("|")+")$" : ''
                }
                api.column(colIdx).search(regex, this.value != '', this.value == '').draw();
            });
    });
}

/**
 * Adds a wrapper around a dataTable to enable horizontal scrolling. This is handled in this
 * function instead of with css because:
 *  - wrapping the actual <table> element will also scroll the dataTables search bar
 *  - the dataTables `enableScrollX` option is buggy and doesn't work
 * @param {String} querySelector a querySelector for this element (ex. '#table_id')
 */
export function enableScrollX(querySelector){
    $(`.dataTable${querySelector}`).wrap("<div style='width: 100%; overflow-x:scroll'></div>");
}

/**
 * This function is used in the exportOptions of a DataTable's "button" configuration as follows:
 * exportOptions: {
 *   format: {
 *     body: (data, row, column, node) => dtExcelExportFormatBody(data, node)
 *   }
 * }
 * This modifies data from table cells that have elements in them (<select>,<a>, etc.) so that
 * the data in the export is formatted properly (ex: an Excel export of a <select> field should have
 * the text from the chosen option)
 * @param {String} data the current cell's HTML content as a string (excluding <td></td>)
 * @param {Object} node the current cell's HTML element (full <td>)
 * @returns 
 */
export function dtExcelExportFormatBody(data, node) {
    const selectElements = node.getElementsByTagName('select')
    const linkElements = node.getElementsByTagName('a')
    const textAreaElements = node.getElementsByTagName("textarea")
    const divElements = node.getElementsByTagName("div")
    if(selectElements.length > 0){
        return selectElements[0].options[selectElements[0].selectedIndex].text
    }else if(linkElements.length > 0){
        return linkElements[0].innerText
    }else if(textAreaElements.length > 0){
        return textAreaElements[0].value
    }else if(divElements.length > 0){
        return divElements[0].innerText
    }else{
        return data
    }
}
