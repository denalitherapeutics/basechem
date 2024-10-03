/**
 * Displays the loading modal. If a `formId` is passed in, only display the modal if 
 * all the required fields of the form have values. This is because missing required fields
 * do not trigger a page refresh, so if there are required fields missing values, the loading
 * modal will display (disabling the rest of the screen) despite not having submitted the form.
 * @param {String} title text to display above the loading spinner
 * @param {String} message text to display below the loading spinner
 * @param {String} formId the id of the form tag whose submission is triggering this modal
 */
export function showLoadingModal({ title = null, message = null, formId = null } = {}) {
    window.addEventListener('pageshow', (event) => {
        // handles b/f events to not show the loading modal if the page is simply restored
        if (event.persisted) {
            hideLoadingModal()
        }
    })

    if (title) {
        $("#loading-modal-title").text(title)
    }
    if (message) {
        $("#loading-modal-message").text(message)
    }
    if (!formId) {
        $('#loading-modal').modal('show');
    }
    else if (document.getElementById(formId).checkValidity()) {
        // only show the loading modal if all form fields have valid values
        $('#loading-modal').modal('show');
    }
}

/**
 * Hides the loading modal.
 */
export function hideLoadingModal() {
    $('#loading-modal').modal('hide');
}
