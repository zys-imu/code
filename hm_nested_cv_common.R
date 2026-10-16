# Shared functions for nested CV of histone-modification changes.

options(stringsAsFactors = FALSE)

load_hm_packages <- function() {
  needed <- c("randomForest", "e1071", "pROC")
  missing <- needed[!vapply(needed, requireNamespace, logical(1), quietly = TRUE)]
  if (length(missing)) stop("Install required packages: ", paste(missing, collapse = ", "), call. = FALSE)
}

check_file <- function(path, label = "File") {
  if (!file.exists(path)) stop(label, " not found: ", path, call. = FALSE)
}

read_gene_label <- function(path, label) {
  check_file(path, label)
  x <- read.table(path, header = TRUE, sep = "\t", quote = "", comment.char = "",
                  check.names = FALSE, stringsAsFactors = FALSE)
  if (!nrow(x)) stop(label, " is empty: ", path, call. = FALSE)
  id_col <- match("geneName", names(x)); if (is.na(id_col)) id_col <- 1L
  ids <- as.character(x[[id_col]])
  ids <- ids[!is.na(ids) & nzchar(ids)]
  unique(ids)
}

read_signal <- function(path, sample_name) {
  check_file(path, paste0("Signal file for ", sample_name))
  x <- read.table(path, header = TRUE, sep = "\t", quote = "",
                  comment.char = "", check.names = FALSE, stringsAsFactors = FALSE)
  if (ncol(x) < 2L) stop("Signal file has fewer than 2 columns: ", path, call. = FALSE)
  id_col <- match("geneName", names(x)); if (is.na(id_col)) id_col <- 1L
  ids <- as.character(x[[id_col]])
  value_cols <- setdiff(seq_len(ncol(x)), id_col)
  value_cols <- value_cols[vapply(x[value_cols], function(z) {
    all(is.finite(suppressWarnings(as.numeric(z))))
  }, logical(1))]
  if (!length(value_cols)) stop("No numeric signal columns in: ", path, call. = FALSE)
  names(x)[id_col] <- "gene_id"
  out <- x[, c(id_col, value_cols), drop = FALSE]
  out$gene_id <- ids
  if (anyDuplicated(out$gene_id)) {
    out <- aggregate(. ~ gene_id, out, mean)
  }
  out
}

sample_path <- function(signal_dir, cell, hm, replicate) {
  # Expected filenames: <cell>_<HM>_<replicate>, optionally with .bed.
  base <- file.path(signal_dir, sprintf("%s_%s_%d", cell, hm, replicate))
  candidates <- c(base, paste0(base, ".bed"), paste0(base, ".txt"))
  hit <- candidates[file.exists(candidates)]
  if (!length(hit)) stop("No signal file found; tried: ", paste(candidates, collapse = ", "), call. = FALSE)
  hit[[1L]]
}

average_cell_signal <- function(signal_dir, cell, hm, replicates = 1:2) {
  tabs <- lapply(replicates, function(r) read_signal(sample_path(signal_dir, cell, hm, r), paste(cell, hm, r)))
  ids <- Reduce(intersect, lapply(tabs, `[[`, "gene_id"))
  if (!length(ids)) stop("No common genes across replicates for ", cell, "/", hm, call. = FALSE)
  mats <- lapply(tabs, function(z) {
    z <- z[match(ids, z$gene_id), , drop = FALSE]
    as.matrix(z[, setdiff(names(z), "gene_id"), drop = FALSE])
  })
  if (length(unique(vapply(mats, ncol, integer(1)))) != 1L)
    stop("Replicate signal files have different numbers of windows: ", cell, "/", hm, call. = FALSE)
  avg <- Reduce(`+`, mats) / length(mats)
  colnames(avg) <- colnames(mats[[1L]])
  data.frame(gene_id = ids, avg, check.names = FALSE)
}

make_hm_changes <- function(signal_dir, hms, cells = c("hepatocyte", "HepG2"), replicates = 1:2,
                             pseudocount = 0.01) {
  result <- list()
  for (hm in hms) {
    ref <- average_cell_signal(signal_dir, cells[[1L]], hm, replicates)
    case <- average_cell_signal(signal_dir, cells[[2L]], hm, replicates)
    ids <- intersect(ref$gene_id, case$gene_id)
    ref <- ref[match(ids, ref$gene_id), , drop = FALSE]
    case <- case[match(ids, case$gene_id), , drop = FALSE]
    refv <- as.matrix(ref[, -1L, drop = FALSE]); casev <- as.matrix(case[, -1L, drop = FALSE])
    if (!identical(colnames(refv), colnames(casev))) stop("Window columns differ for HM: ", hm, call. = FALSE)
    z <- log2((casev + pseudocount) / (refv + pseudocount))
    colnames(z) <- paste0(hm, "_", colnames(z))
    result[[hm]] <- data.frame(gene_id = ids, z, check.names = FALSE)
  }
  result
}

make_labels <- function(up_file, down_file) {
  up <- read_gene_label(up_file, "Up-DEG file"); down <- read_gene_label(down_file, "Down-DEG file")
  if (length(intersect(up, down))) stop("Genes occur in both Up and Down DEG lists.", call. = FALSE)
  data.frame(gene_id = c(up, down), y = factor(c(rep("Up", length(up)), rep("Down", length(down))), levels = c("Down", "Up")))
}

stratified_folds <- function(y, k = 5L, seed = 1L) {
  set.seed(seed); folds <- vector("list", k)
  for (lev in levels(y)) {
    idx <- which(y == lev); idx <- sample(idx)
    for (i in seq_along(idx)) folds[[((i - 1L) %% k) + 1L]] <- c(folds[[((i - 1L) %% k) + 1L]], idx[[i]])
  }
  lapply(folds, sort)
}

fit_predict <- function(train_x, train_y, test_x, method) {
  if (method == "RF") {
    fit <- randomForest::randomForest(x = train_x, y = train_y)
    return(as.numeric(predict(fit, test_x, type = "prob")[, "Up"]))
  }
  if (method == "SVM") {
    fit <- e1071::svm(x = train_x, y = train_y, probability = TRUE)
    pr <- attr(predict(fit, test_x, probability = TRUE), "probabilities")
    return(as.numeric(pr[, "Up"]))
  }
  fit <- stats::glm(y ~ ., data = data.frame(y = train_y, train_x), family = stats::binomial())
  as.numeric(stats::predict(fit, newdata = data.frame(test_x), type = "response"))
}

auc_score <- function(y, p) {
  if (length(unique(y)) < 2L || any(!is.finite(p))) return(NA_real_)
  as.numeric(pROC::auc(pROC::roc(response = y, predictor = p, levels = c("Down", "Up"), direction = "<", quiet = TRUE)))
}

classification_metrics <- function(y, p, threshold = 0.5) {
  pred <- factor(ifelse(p >= threshold, "Up", "Down"), levels = c("Down", "Up"))
  tp <- sum(pred == "Up" & y == "Up"); tn <- sum(pred == "Down" & y == "Down")
  fp <- sum(pred == "Up" & y == "Down"); fn <- sum(pred == "Down" & y == "Up")
  data.frame(AUC = auc_score(y, p), sensitivity = if ((tp + fn) > 0) tp / (tp + fn) else NA_real_,
             specificity = if ((tn + fp) > 0) tn / (tn + fp) else NA_real_)
}

inner_auc <- function(x, y, features, method, folds) {
  scores <- vapply(folds, function(test_idx) {
    train_idx <- setdiff(seq_along(y), test_idx)
    p <- tryCatch(fit_predict(x[train_idx, features, drop = FALSE], y[train_idx], x[test_idx, features, drop = FALSE], method), error = function(e) rep(NA_real_, length(test_idx)))
    auc_score(y[test_idx], p)
  }, numeric(1))
  mean(scores, na.rm = TRUE)
}

select_greedy <- function(x, y, method, folds, max_features = 10L) {
  remaining <- colnames(x); selected <- character(); best <- -Inf
  while (length(remaining) && length(selected) < max_features) {
    scores <- vapply(remaining, function(f) inner_auc(x, y, c(selected, f), method, folds), numeric(1))
    if (!any(is.finite(scores))) break
    j <- which.max(scores); if (scores[[j]] <= best + 1e-8) break
    selected <- c(selected, remaining[[j]]); remaining <- remaining[-j]; best <- scores[[j]]
  }
  if (!length(selected)) selected <- colnames(x)[1L]
  selected
}

select_subset <- function(x, y, method, folds, max_subset_size = ncol(x)) {
  features <- colnames(x); max_subset_size <- min(max_subset_size, length(features)); best <- -Inf; winner <- features[1L]
  for (size in seq_len(max_subset_size)) {
    combos <- combn(features, size, simplify = FALSE)
    for (set in combos) { score <- inner_auc(x, y, set, method, folds); if (is.finite(score) && score > best) { best <- score; winner <- set } }
  }
  winner
}

nested_oof <- function(x, y, method, selector = c("greedy", "subset", "fixed"), outer_k = 10L, inner_k = 5L, seed = 1L, max_features = 10L, max_subset_size = ncol(x)) {
  selector <- match.arg(selector); outer <- stratified_folds(y, outer_k, seed); p <- rep(NA_real_, length(y)); selected <- vector("list", length(outer))
  for (i in seq_along(outer)) {
    test_idx <- outer[[i]]; train_idx <- setdiff(seq_along(y), test_idx); inner <- stratified_folds(y[train_idx], inner_k, seed + i)
    inner <- lapply(inner, function(z) train_idx[z]);
    feats <- if (selector == "fixed") colnames(x) else if (selector == "greedy") select_greedy(x[train_idx, , drop = FALSE], y[train_idx], method, inner, max_features) else select_subset(x[train_idx, , drop = FALSE], y[train_idx], method, inner, max_subset_size)
    selected[[i]] <- feats
    p[test_idx] <- tryCatch(fit_predict(x[train_idx, feats, drop = FALSE], y[train_idx], x[test_idx, feats, drop = FALSE], method), error = function(e) NA_real_)
  }
  list(prediction = p, metrics = classification_metrics(y, p), selected = selected)
}
