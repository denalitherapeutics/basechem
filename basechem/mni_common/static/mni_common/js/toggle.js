/**
 * Disables/Enables all child elements within the given element
 * @param {Object} element an html element
 * @param {Boolean} disable true if the element should be disabled, false if it should be enabled
 * @param {Boolean} buttons true if this should also disable/enable buttons, false if buttons should be ignored
 */
export function toggleDisableChildren(element, disable, buttons=true){
    const children = element.getElementsByTagName('*');
    for(let i = 0; i < children.length; i++){
        if(buttons || children[i].tagName != "BUTTON" ){
          children[i].disabled = disable;
        }
    }
}

/**
 * Displays or hides an html element. If the element has a label and formInput (i.e. it is
 * a form element), updates the "required" status of the field based on the value of `required`.
 * @param {String} element_id the id of an html element
 * @param {Boolean} show true to display, false to hide
 * @param {Boolean} required a boolean, should this form field be required. This param is only
 *  relevant for form fields when `show` is true. If the element is not a form field or `show`
 *  is false, the element will not be required even if this param is true.
 */
 export function toggleElement(element_id, show, required=false){
    // Show/hide element
    const element = document.getElementById(element_id)
    element.style.display = show ? '' : 'none';
  
    const label = element.querySelector('label');
    var formInput = element.querySelector('.form-control')
    if(!formInput){
      // Checkbox fields have a different control class
      formInput = element.querySelector('.custom-control-input')
    }
    
    // Show/hide required label if element is a form field
    if(label && formInput){
      if(show && required){
        label.classList.add("requiredField")
        // Add asterisk field to label
        if(!label.querySelector(".asteriskField")){
          const asteriskField = document.createElement("span");
          asteriskField.innerText = "*"
          asteriskField.classList.add("asteriskField")
          label.appendChild(asteriskField);
        }
        // Add `required` to input field
        formInput.required = true;
      }else{
        label.classList.remove("requiredField")
        // Remove asterisk field if it exists
        const asteriskField = label.querySelector('.asteriskField')
        if(asteriskField){
          label.removeChild(asteriskField);
        }
        // Remove `required` from input field
        formInput.required = false;
  
      }
    }
  }
