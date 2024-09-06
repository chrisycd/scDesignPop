################# Helper functions #################

#' Check if multiple vectors have same elements
#'
#' This is an internal helper function to check if two or more vectors have same
#'     elements (with their order considered as an option)
#'
#' @param ... arguments of two or more vectors
#' @param ignore_order logical scalar to disregard the ordering of elements in
#'     input vectors. Default is TRUE.
#'
#' @return a logical scalar
#' @export
#'
#' @examples
#' vec1 <- c("cherry", "apple", "banana", "cherry")
#' vec2 <- c("cherry", "apple", "cherry", "banana")
#' vec3 <- c("banana", "cherry", "apple", "cherry")
#'
#' checkVectorEqual(vec1, vec2, vec3, ignore_order = TRUE)  # returns TRUE
#' checkVectorEqual(vec1, vec2, vec3, ignore_order = FALSE)  # returns FALSE
checkVectorEqual <- function(..., ignore_order = TRUE) {

    vec_list <- list(...)

    if(length(vec_list) < 2) {
        stop(sprintf("Input must be 2 or more vectors!"))
    }

    first_vec <- vec_list[[1]]

    # compare all vectors
    if(ignore_order) {
        res <- all(base::sapply(vec_list[-1], function(x) base::setequal(first_vec, x)))
    } else {
        res <- all(base::sapply(vec_list[-1], function(x) base::identical(first_vec, x)))
    }

    return(res)
}


#' Check membership of first vector compared to other vectors
#'
#' @param ... arguments of two or more vectors
#'
#' @return a logical scalar
#' @export
#'
#' @examples
#' NULL
checkVectorContain <- function(...) {

    vec_list <- list(...)

    if(length(vec_list) < 2) {
        stop(sprintf("Input must be 2 or more vectors!"))
    }

    first_vec <- vec_list[[1]]

    # compare membership of first vector to all others
    res <- all(base::sapply(vec_list[-1], function(x) base::is.element(first_vec, x)))

    return(res)
}
