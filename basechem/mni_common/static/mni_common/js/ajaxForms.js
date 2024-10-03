import { hideLoadingModal } from "./loading.js"

/**
 * Handles form submission for all forms with the class "ajax-form". Any `form` element
 * with this class needs to have the following specified:
 * - form-errors: the id of an empty div element where form-level errors should be displayed
 * - action: a url to post to (usually the django URL of the form view)
 * - method: usually "POST"
 * 
 * An example: <form class="ajax-form" id="ex-form" form-errors="ex-form-errors" method="POST" action="{% url ex_form_view my_details %}"
 * The form view should return a json response
 */
export default function(){
  $('.ajax-form').on('submit', function(event){
    event.preventDefault();
    postAjaxForm.call(this);
  });
}

/**
 * Called on a form element of the class `ajax-form`, this handles the form post, displays
 * error and success messages, and redirects the page as needed. The expected response from the
 * post request is a JSON object with the following fields
 * - {Object} errors: a dictionary of errors to display to the user of the form {"__all__": [form level errors], "field1": [field1 errors], etc.}
 * - {String} success_method: optional, a javascript method to run if there are no errors (for example, to update other elements on the page)
 * - {String} redirect_url: optional, the url to redirect to if there are no errors
 * - {Boolean} close_modal: optional, true if the modal containing the form should close if there are no errors (instead of reloading)
 * 
 * If redirect_url or close_modal are specified, the page will redirect/close the modal instead of reloading. Otherwise, the default behavior
 * is to reload the page.
 */
export function postAjaxForm() {
  const errorDiv = document.getElementById($(this).attr('form-errors'))
  const formData = new FormData(this);
  const form = $(this)[0]
  toggleForm(form, false)
  const formLoadingSpinner = form.querySelector("#form-loading-spinner")
  // If there is a loading spinner inside the form element, display it now
  if(formLoadingSpinner){
    formLoadingSpinner.style.display = ""
  }
  // If this form has a prefix (to avoid clashes with other forms on the page), remove
  // the prefix from the form data, as django expects to see the field names without prefixes
  if(form.hasAttribute("form-prefix")){
    const formKeys = Array.from(formData.keys())
    const prefix = form.getAttribute("form-prefix")
    for(const key of formKeys){
      if(key.includes(prefix)){
        const newKey = key.replace(prefix+"-", "")
        for(const value of formData.getAll(key)){
          formData.append(newKey, value)
        }
        formData.delete(key)
      }
    }
  }
  $.ajax({
    url : this.action,
    type : $(this).attr('method'),
    data : formData,
    processData: false,
    contentType: false,
    success : function(response) {
      // If there is a loading spinner inside the form element, hide it 
      if(formLoadingSpinner){
        formLoadingSpinner.style.display = "none"
      }
      if(!Object.keys(response.errors).length){
        // No errors, handle form success
        clearFormMessages(form)
        if(response.success_method){
          eval(response.success_method)
        }
        if(response.redirect_url) {
          window.location = response.redirect_url
        }else if(response.close_modal){
          closeModal(form)
          toggleForm(form,true)
          // If a loading modal was displayed, hide it
          hideLoadingModal()
        }else{
          window.location.reload()
        }
      }else{
        displayErrors(response.errors, errorDiv, form)
        toggleForm(form,true)
        // If a loading modal was displayed, hide it
        hideLoadingModal()
      }
    },error: function(){
      const message = `<div class="alert alert-block alert-danger">An unexpected error occurred. Please contact MnI</div>`
      errorDiv.innerHTML += message
      // Scroll the modal so the general errors are visible
      errorDiv.scrollIntoView({ behavior: 'smooth', block: 'nearest', inline: 'start' })
      // If there is a loading spinner inside the form element, hide it 
      if(formLoadingSpinner){
        formLoadingSpinner.style.display = "none"
      }
    }
  });
};
  
/**
 * Given a json object returned from an ajax form post, display error messages to the user. This
 * function makes some assumptions about the structure of a django form (specifically that
 * each input field's id is "id_field_name"). If you need to override field names (for example,
 * when there are multiple copies of the form on a single page), the field id `formId_fieldName`
 * also works.
 * @param {Object} errors an object of the form {field1: error1, field2: error2, "__all__": [errors]}
 *  where __all__ contains a list of errors that are not associated with a particular field
 * @param {Object} errorDiv an html element where form-level errors should be displayed
 * @param {String} form the form element whose errors are being displayed
 */
function displayErrors(errors, errorDiv, form){
  clearFormMessages(form)
  for(const fieldKey in errors){
    if(fieldKey!="__all__"){
      addFieldError(form, fieldKey, errors[fieldKey])
    }else{
      addFormError(form, errorDiv, errors.__all__)
    }
  }
}

/**
 * Use this function with a `hidden.bs.modal` event trigger to reset the form when
 * a modal is closed (to prevent a user from thinking they've saved information that they haven't).
 * 
 * Example call:
 * $(modal-id).on("hidden.bs.modal", resetForm)
 * @param {Object} ev event object 
 */
export const resetForm = ev => {
  const modal = $(ev.currentTarget)[0]
  const form = $(modal).find("form")[0]
  form.reset()
  clearFormMessages(form)
}

/*************************
 *         UTILS         *
 *************************/

/**
 * Clears all messages on a form including form-level error or success messages and field-level
 * error messages
 * @param {Object} form a form element
 */
export function clearFormMessages(form){
  form.querySelectorAll('.invalid-feedback').forEach(e => e.remove());
  form.querySelectorAll('.is-invalid').forEach(e => e.classList.remove("is-invalid"))
  form.querySelectorAll('.alert.alert-block.alert-danger').forEach(e => e.remove());
  form.querySelectorAll('.alert.alert-block.alert-success').forEach(e => e.remove());
}

/**
 * Enables/Disables all the input fields and buttons in a form
 * @param {Object} form the form that should be enabled/disabled
 * @param {Boolean} enable true if the form should be enabled, false if it should be disabled
 */
function toggleForm(form, enable){
  if(enable){
    // Enable form
    $(form).find('*').removeAttr('disabled');
    form.querySelectorAll("iframe").forEach( iframe => iframe.style.pointerEvents = "auto")
  }else{
    // Disable form
    $(form).find('*').attr('disabled', true);
    form.querySelectorAll("iframe").forEach( iframe => iframe.style.pointerEvents = "none")
  }
}

/**
 * Closes the parent modal of form.
 * @param {String} form the form element whose parent modal should be closed
 */
function closeModal(form){
  const modal = form.closest(".modal")
  $(`#${modal.id}`).modal('hide');
}

/**
 * Clears the currently selected file from the given file input. This is needed because of a 
 * chromium #wontfix bug (1086707). If a user (1)chooses a file, (2)submits the form & gets back errors,
 * (3)edits the same file w/ same filename, and (4) resubmits the form, chromium raises
 * an exception because of a potential security risk (there's a lot of debate about this).
 * 
 * This function avoids that by clearing the file input so the user is required to re-upload
 * the file and there is no cached version of the file stored in the browser.
 * @param {Object} form the form element whose file field needs to be cleared
 * @param {Object} field the file input element
 */
function clearFileField(form, field) {
  field.value = ''
  // The label is used by browsers to search for the file, so this also needs to be cleared
  const label = form.querySelector(`.custom-file-label[for=${field.id}]`)
  label.innerText = ''
}

/**
 * Given a form and a field name (ex: "project_name"), return the element for this input
 * field in the given form
 * @param {Object} form a form element
 * @param {String} fieldKey the name of the field
 * @returns a field element in the given form (if found)
 */
function getFieldElement(form, fieldKey){
  if(form.hasAttribute("form-prefix")){
    var fieldElement = form.querySelector(`#id_${form.getAttribute("form-prefix")}-${fieldKey}`)
  }else{
    var fieldElement = form.querySelector(`#id_${fieldKey}`)
  }
  if(!fieldElement){
    // For forms that override the default IDs, the field id might be
    // prepended by the form ID instead of 'id'
    fieldElement = form.querySelector(`#${form.id}_${fieldKey}`)
  }
  return fieldElement
}

/**
 * Given a form, field name, and error messages, add the given error messages to the field
 * @param {Object} form a form element
 * @param {String} fieldKey the name of the field whose errors are passed
 * @param {Array} errors an array of error messages (strings) for this field
 */
function addFieldError(form, fieldKey, errors){
  const fieldElement = getFieldElement(form, fieldKey)
  // Add red border to field
  fieldElement.classList.add("is-invalid")
  for(var i=0; i<errors.length; i++){
    // Create container for error messages
    const errorContainer = document.createElement('span');
    errorContainer.id=`error_${i}_id_${fieldKey}`
    errorContainer.classList.add("invalid-feedback")
    // Add error message to container
    var errorMessage = document.createElement('strong');
    errorMessage.innerText = errors[i]
    // Append error messages below input field
    errorContainer.appendChild(errorMessage)
    fieldElement.after(errorContainer)
    // Clear field if file input
    if (fieldElement.getAttribute("type") == "file") {
      clearFileField(form, fieldElement)
    }
  }
}

/**
 * Add the given error messages to the form and scroll the errorDiv into view
 * @param {Object} form a form element
 * @param {Object} errorDiv a div element that should contain the errors for this form
 * @param {Array} errors a list of form-level error messages
 */
function addFormError(form, errorDiv, errors){
    var message = `<div class="alert alert-block alert-danger">Please correct the following errors:<ul style="margin:0px">`
    for(var i = 0; i < errors.length; i++){
      message = message + `<li>${errors[i]}</li>`
    }
    errorDiv.innerHTML += message +'</ul></div>'
    form.querySelectorAll('input[type=file]').forEach(field => clearFileField(form, field));
    // Scroll the modal so the form-levels errors are in-view
    errorDiv.scrollIntoView({ behavior: 'smooth', block: 'nearest', inline: 'start' })
}