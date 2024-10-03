/**
 * Copies an element to the user's clipboard
 * @param {String} querySelector a unique querySelector for the element that is being copied
 * @param {Boolean} asImage should this element be copied as a png? If False, the element's innerText is copied instead
 */
export function copyToClipboard(querySelector, asImage) {
    // Check if the user is using Chrome because not all browsers allow writing images to the clipboard
    const isChrome = navigator.userAgent.indexOf("Chrome") > -1;
    const element = document.querySelector(querySelector);
    if (asImage && isChrome) {
        // Convert element to png and copy image to clipboard
        htmlToImage.toBlob(element, { skipFonts: true }).then(async function (blob) {
            navigator.clipboard.write([new ClipboardItem({ [blob.type]: blob })])
        })
    } else {
        // Copy the text inside the element to the clipboard
        navigator.clipboard.writeText(element.innerText);
    }
}
