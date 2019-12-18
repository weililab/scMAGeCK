# function definitions ##### version: 02-22-2019 should sync from the version in macbook:
# /Users/weili/Dropbox/work/cropseq/Shendure/nmeth18/multiple_guides_function.R


getscaledata <- function(targetobj, scaled = TRUE) {
    # if scaled=FALSE, return raw.data
    if ("scale.data" %in% names(attributes(targetobj))) {
        if (scaled) {
            scalef = targetobj@scale.data  # for version 2
        } else {
            scalef = targetobj@raw.data  # for version 2
        }
    } else {
        if (scaled) {
            scalef = GetAssayData(object = targetobj, slot = "scale.data")
        } else {
            scalef = GetAssayData(object = targetobj, slot = "counts")
        }
    }
    return(scalef)
}
TRUE
