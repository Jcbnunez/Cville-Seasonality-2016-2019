sumz <-
function (p, weights = NULL, data = NULL, subset = NULL, na.action = na.fail,
    log.p = FALSE)
{
    if (is.null(data))
        data <- sys.frame(sys.parent())
    mf <- match.call()
    mf$data <- NULL

    mf$subset <- NULL
    mf$na.action <- NULL
    mf[[1]] <- as.name("data.frame")
    mf <- eval(mf, data)
    if (!is.null(subset))
        mf <- mf[subset, ]
    mf <- na.action(mf)
    p <- as.numeric(mf$p)
    weights <- mf$weights
    noweights <- is.null(weights)
    if (noweights)
        weights <- rep(1, length(p))
    if (length(p) != length(weights))
        warning("Length of p and weights differ")
    keep <- (p > 0) & (p < 1)
    invalid <- sum(1L * keep) < 2
    if (invalid) {
        warning("Must have at least two valid p values")
        res <- list(z = NA_real_, p = NA_real_, validp = p[keep],
            weights = weights)
    }

    else {
        if (sum(1L * keep) != length(p)) {
            warning("Some studies omitted")
            omitw <- weights[!keep]
            if ((sum(1L * omitw) > 0) & !noweights)
                warning("Weights omitted too")
        }
        zp <- (qnorm(p[keep], lower.tail = FALSE) %*% weights[keep])/sqrt(sum(weights[keep]^2))
        res <- list(z = zp, p = pnorm(zp, lower.tail = FALSE,
            log.p = log.p), validp = p[keep], weights = weights)
    }
    class(res) <- c("sumz", "metap")
    res$z
}
