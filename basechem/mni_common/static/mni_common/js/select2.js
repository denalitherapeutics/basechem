/**
 * Used in select2 configs, this copies the `data-id` attribute from the html select element
 * to the select2 element that replaces it. This allows us to use `dataId` as an identifier
 * for javascript that interacts with select2 elements. This is needed because custom html
 * attributes don't get copied to the select2 html that replaces the original html.
 * 
 * EXAMPLE:
 * html: `<select id="display-reference"><option data-id="s-1" value="s1"/></select>`
 *
 * select2 config: `$('#display-reference').select2({
      templateResult: copyDataId,
      templateSelection: copyDataId
    })`

 * resulting html: `<span class="select2-selection__rendered" id="select2-display-receptor-container" data-id="s-1"></span>`
 * without this, the resulting html would be missing `data-id="s-1"`
 */
export function copyDataId (data, container){
    if(data.element){
        if(data.element.getAttribute("data-id")){
        $(container).attr("data-id", data.element.getAttribute("data-id"));
        }
        // If there is display text, use that instead of the default text
        const display = $(data.element).data("display")
        if(display){
            return display
        }
    }
    return data.text
}
