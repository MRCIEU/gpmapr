.prepare_snp_similarity <- function(s_matrix, similarity_threshold = NULL) {
  if (!is.matrix(s_matrix) || nrow(s_matrix) != ncol(s_matrix)) {
    stop("s_matrix must be a square matrix")
  }
  if (!is.null(similarity_threshold) && similarity_threshold < 0) {
    stop("similarity_threshold must be NULL, 0, or a non-negative value")
  }

  snp_ids <- colnames(s_matrix)
  if (is.null(snp_ids)) {
    snp_ids <- as.character(seq_len(ncol(s_matrix)))
  }

  S <- s_matrix
  dimnames(S) <- list(snp_ids, snp_ids)
  diag(S) <- 1
  S[!is.finite(S)] <- 0

  if (!is.null(similarity_threshold) && similarity_threshold > 0) {
    S <- .sparsify_similarity_matrix(S, threshold = similarity_threshold)
  }

  return(list(s_matrix = S, snp_ids = snp_ids))
}


.sparsify_similarity_matrix <- function(s_matrix, threshold) {
  n <- nrow(s_matrix)
  out <- s_matrix
  off_diag <- matrix(TRUE, n, n)
  diag(off_diag) <- FALSE
  out[off_diag & abs(out) < threshold] <- 0
  diag(out) <- 1
  return(out)
}


#' @title Cluster SNPs With Signed Louvain
#' @description Cluster SNPs from a cosine similarity matrix by maximising Gomez
#' signed modularity. Negative edges (anti-parallel profiles) are retained.
#' The number of modules is data-driven for fixed `gamma` and
#' `similarity_threshold`.
#' @param s_matrix Symmetric SNP-by-SNP cosine similarity matrix from
#'   `snp_similarity_matrix()$s_matrix`.
#' @param similarity_threshold Minimum absolute similarity for off-diagonal SNP
#'   pairs. Weaker edges are zeroed before clustering. Set to `NULL` or `0` to
#'   use the full matrix. Defaults to `NULL`.
#' @param gamma Resolution parameter. Values greater than 1 tend to yield more,
#'   smaller modules. Defaults to `1`.
#' @param seed Optional RNG seed for Louvain node-order randomisation.
#'   Defaults to `1L`.
#' @return A list with:
#'   \itemize{
#'     \item cluster: named integer vector of cluster assignments per SNP
#'     \item method: `"louvain"`
#'     \item n_clusters: number of clusters returned
#'     \item similarity_threshold: edge threshold applied before clustering
#'     \item gamma: Louvain resolution used
#'     \item details: signed graph and modularity diagnostics
#'   }
#' @export
cluster_snp_profiles_louvain <- function(s_matrix,
                                         similarity_threshold = NULL,
                                         gamma = 1,
                                         seed = 1L) {
  if (gamma <= 0) {
    stop("gamma must be positive")
  }

  prepared <- .prepare_snp_similarity(
    s_matrix,
    similarity_threshold = similarity_threshold
  )
  S <- prepared$s_matrix
  snp_ids <- prepared$snp_ids

  affinity <- S
  diag(affinity) <- 0
  affinity[!is.finite(affinity)] <- 0

  signed <- .cluster_louvain_signed(affinity, gamma = gamma, seed = seed)
  graph <- .igraph_from_signed_adjacency(affinity)
  clusters <- as.integer(signed$cluster)
  names(clusters) <- snp_ids

  return(list(
    cluster = clusters,
    method = "louvain",
    n_clusters = length(unique(clusters)),
    similarity_threshold = similarity_threshold,
    gamma = gamma,
    details = list(
      igraph = graph,
      modularity = signed$modularity,
      qtype = signed$qtype,
      gamma = gamma,
      algorithm = "signed_louvain"
    )
  ))
}


#' @title Cluster SNPs With Spectral Clustering
#' @description Cluster SNPs from a cosine similarity matrix using the normalised
#' graph Laplacian of the shifted affinity matrix \eqn{(S + 1) / 2}, followed by
#' k-means in the eigenvector embedding.
#' @param s_matrix Symmetric SNP-by-SNP cosine similarity matrix from
#'   `snp_similarity_matrix()$s_matrix`.
#' @param k Number of clusters. Must be at least 2.
#' @param similarity_threshold Minimum absolute similarity for off-diagonal SNP
#'   pairs. Weaker edges are zeroed before clustering. Set to `NULL` or `0` to
#'   use the full matrix. Defaults to `NULL`.
#' @return A list with:
#'   \itemize{
#'     \item cluster: named integer vector of cluster assignments per SNP
#'     \item method: `"spectral"`
#'     \item k: requested number of clusters
#'     \item n_clusters: number of clusters returned
#'     \item similarity_threshold: edge threshold applied before clustering
#'     \item details: embedding, k-means fit, and leading eigenvalues
#'   }
#' @export
cluster_snp_profiles_spectral <- function(s_matrix,
                                          k = 3,
                                          similarity_threshold = NULL) {
  if (k < 2) {
    stop("k must be at least 2")
  }

  prepared <- .prepare_snp_similarity(
    s_matrix,
    similarity_threshold = similarity_threshold
  )
  S <- prepared$s_matrix
  snp_ids <- prepared$snp_ids

  affinity <- (S + 1) / 2
  affinity[affinity < 0] <- 0
  diag(affinity) <- 0

  degree <- rowSums(affinity)
  degree[degree == 0] <- 1
  d_inv_sqrt <- diag(1 / sqrt(degree))
  laplacian <- diag(nrow(affinity)) - d_inv_sqrt %*% affinity %*% d_inv_sqrt

  eig <- eigen(laplacian, symmetric = TRUE)
  embedding <- eig$vectors[, seq_len(k), drop = FALSE]
  row_norms <- sqrt(rowSums(embedding^2))
  row_norms[row_norms == 0] <- 1
  embedding <- embedding / row_norms

  km <- stats::kmeans(embedding, centers = k, nstart = 25)
  clusters <- km$cluster
  names(clusters) <- snp_ids

  return(list(
    cluster = clusters,
    method = "spectral",
    k = k,
    n_clusters = length(unique(clusters)),
    similarity_threshold = similarity_threshold,
    details = list(
      embedding = embedding,
      kmeans = km,
      eigenvalues = eig$values[seq_len(k)]
    )
  ))
}


.igraph_from_signed_adjacency <- function(affinity) {
  adj <- affinity
  diag(adj) <- 0
  adj[!is.finite(adj)] <- 0

  presence <- abs(adj) > 0
  graph <- igraph::graph_from_adjacency_matrix(
    presence * 1,
    mode = "undirected",
    weighted = FALSE,
    diag = FALSE
  )

  if (igraph::ecount(graph) > 0) {
    edge_ends <- igraph::as_edgelist(graph)
    edge_weights <- adj[cbind(edge_ends[, 1], edge_ends[, 2])]
    igraph::E(graph)$weight <- edge_weights
  }

  return(graph)
}


.cluster_louvain_signed <- function(s_matrix, gamma = 1, qtype = "sta", seed = NULL) {
  W <- s_matrix
  diag(W) <- 0
  W[!is.finite(W)] <- 0

  n <- nrow(W)
  W0 <- W * (W > 0)
  W1 <- -W * (W < 0)
  s0 <- sum(W0)
  s1 <- sum(W1)

  qtype <- match.arg(qtype, c("sta", "pos", "smp", "gja", "neg"))

  if (qtype == "smp") {
    d0 <- if (s0 > 0) 1 / s0 else 0
    d1 <- if (s1 > 0) 1 / s1 else 0
  } else if (qtype == "gja") {
    denom <- s0 + s1
    d0 <- if (denom > 0) 1 / denom else 0
    d1 <- d0
  } else if (qtype == "sta") {
    d0 <- if (s0 > 0) 1 / s0 else 0
    d1 <- if ((s0 + s1) > 0) 1 / (s0 + s1) else 0
  } else if (qtype == "pos") {
    d0 <- if (s0 > 0) 1 / s0 else 0
    d1 <- 0
  } else if (qtype == "neg") {
    d0 <- 0
    d1 <- if (s1 > 0) 1 / s1 else 0
  }

  if (s0 == 0) {
    s0 <- 1
    d0 <- 0
  }
  if (s1 == 0) {
    s1 <- 1
    d1 <- 0
  }

  if (!is.null(seed)) {
    set.seed(seed)
  }

  h <- 1L
  nh <- n
  ci <- list(NULL, seq_len(n))
  q <- c(-1, 0)

  while (q[h + 1L] - q[h] > 1e-10) {
    if (h > 300L) {
      stop("Signed Louvain exceeded maximum hierarchy depth", call. = FALSE)
    }

    kn0 <- colSums(W0)
    kn1 <- colSums(W1)
    km0 <- kn0
    km1 <- kn1
    knm0 <- W0
    knm1 <- W1

    m <- seq_len(nh)
    flag <- TRUE
    it <- 0L

    while (flag) {
      it <- it + 1L
      if (it > 1000L) {
        stop("Signed Louvain iteration limit exceeded", call. = FALSE)
      }
      flag <- FALSE

      for (u in sample(nh)) {
        ma <- m[u]
        dQ0 <- (knm0[u, ] + W0[u, u] - knm0[u, ma]) -
          gamma * kn0[u] * (km0 + kn0[u] - km0[ma]) / s0
        dQ1 <- (knm1[u, ] + W1[u, u] - knm1[u, ma]) -
          gamma * kn1[u] * (km1 + kn1[u] - km1[ma]) / s1

        dQ <- d0 * dQ0 - d1 * dQ1
        dQ[ma] <- 0

        if (max(dQ) > 1e-10) {
          flag <- TRUE
          mb <- which.max(dQ)

          knm0[, mb] <- knm0[, mb] + W0[, u]
          knm0[, ma] <- knm0[, ma] - W0[, u]
          knm1[, mb] <- knm1[, mb] + W1[, u]
          knm1[, ma] <- knm1[, ma] - W1[, u]
          km0[mb] <- km0[mb] + kn0[u]
          km0[ma] <- km0[ma] - kn0[u]
          km1[mb] <- km1[mb] + kn1[u]
          km1[ma] <- km1[ma] - kn1[u]
          m[u] <- mb
        }
      }
    }

    h <- h + 1L
    ci[[h + 1L]] <- numeric(n)
    m_factor <- as.integer(as.factor(m))

    for (u in seq_len(nh)) {
      ci[[h + 1L]][ci[[h]] == u] <- m_factor[u]
    }

    nh <- max(m_factor)
    wn0 <- matrix(0, nh, nh)
    wn1 <- matrix(0, nh, nh)
    for (u in seq_len(nh)) {
      for (v in u:nh) {
        idx_u <- m_factor == u
        idx_v <- m_factor == v
        val0 <- sum(W0[idx_u, idx_v, drop = FALSE])
        val1 <- sum(W1[idx_u, idx_v, drop = FALSE])
        wn0[u, v] <- val0
        wn0[v, u] <- val0
        wn1[u, v] <- val1
        wn1[v, u] <- val1
      }
    }
    W0 <- wn0
    W1 <- wn1

    q0 <- sum(diag(W0)) - gamma * sum(W0 %*% W0) / s0
    q1 <- sum(diag(W1)) - gamma * sum(W1 %*% W1) / s1
    q <- c(q, d0 * q0 - d1 * q1)
  }

  membership <- as.integer(as.factor(ci[[h + 1L]]))

  return(list(
    cluster = membership,
    modularity = q[length(q)],
    qtype = qtype
  ))
}

#' @title Summarise SNP Module Phenotype Drivers
#' @description For each SNP cluster, compute the mean oriented pleiotropy profile
#' across member SNPs and rank background traits by absolute effect magnitude to
#' identify functional drivers of each biological module.
#' @param x_matrix Oriented pleiotropy matrix (traits x SNPs), typically
#'   `orient_pleiotropy_matrix()$x_matrix`.
#' @param cluster Named integer vector of cluster assignments per SNP column
#'   (as returned by `cluster_snp_profiles_louvain()$cluster` or
#'   `cluster_snp_profiles_spectral()$cluster`).
#' @param trait_info Optional dataframe with `trait_id` and `trait_name` columns
#'   (as returned by `build_pleiotropy_matrix()$trait_info`).
#' @param exclude_trait_id Optional trait ID to exclude from summaries (e.g. the
#'   target trait when interpreting background drivers).
#' @param na_as_zero Treat `NA` entries as zero before averaging. Defaults to `TRUE`.
#' @param top_n If not `NULL`, return only the top `n` traits per cluster ranked by
#'   `abs_mean_z`.
#' @return A dataframe with columns: `cluster`, `trait_id`, `trait_name` (if
#'   available), `mean_z`, `abs_mean_z`, `n_snps`.
#' @export
summarise_snp_modules <- function(x_matrix,
                                  cluster,
                                  trait_info = NULL,
                                  exclude_trait_id = NULL,
                                  na_as_zero = TRUE,
                                  top_n = NULL) {
  if (length(cluster) != ncol(x_matrix)) {
    stop("length(cluster) must match ncol(x_matrix)")
  }

  if (is.null(names(cluster)) && !is.null(colnames(x_matrix))) {
    names(cluster) <- colnames(x_matrix)
  }

  x <- x_matrix
  if (na_as_zero) {
    x[is.na(x)] <- 0
  }

  cluster_ids <- sort(unique(cluster))
  summaries <- lapply(cluster_ids, function(cl) {
    snp_cols <- names(cluster)[cluster == cl]
    module_matrix <- x[, snp_cols, drop = FALSE]
    mean_z <- rowMeans(module_matrix, na.rm = TRUE)

    out <- data.frame(
      cluster = cl,
      trait_id = rownames(x_matrix),
      mean_z = mean_z,
      abs_mean_z = abs(mean_z),
      n_snps = length(snp_cols),
      stringsAsFactors = FALSE
    )

    if (!is.null(exclude_trait_id)) {
      out <- out[out$trait_id != as.character(exclude_trait_id), , drop = FALSE]
    }

    out <- out[order(-out$abs_mean_z), , drop = FALSE]

    if (!is.null(top_n)) {
      out <- utils::head(out, top_n)
    }

    return(out)
  })

  result <- dplyr::bind_rows(summaries)

  if (!is.null(trait_info)) {
    trait_info <- trait_info |>
      dplyr::mutate(trait_id = as.character(trait_id))
    result <- result |>
      dplyr::left_join(trait_info, by = "trait_id")
  }

  return(result)
}


#' @title SNP Module Reliability Metrics
#' @description Graph-quality diagnostics for SNP modules on a cosine similarity
#' matrix. For each module this reports internal vs external similarity,
#' connectedness after applying `edge_threshold` to `|S|`, and a diagnostic
#' silhouette. Use size, internal similarity, and connectedness to drop weakly
#' linked modules rather than keeping every Louvain community.
#'
#' Connectedness is the fraction of unique SNP pairs in the module with
#' \eqn{|S| \ge} `edge_threshold`. Silhouette uses distance \eqn{1 - S} and is
#' reported only as a diagnostic: because SNPs can participate in multiple
#' biological programs, negative silhouettes are expected and are not used to
#' gate `reliable`.
#' @param s_matrix Symmetric SNP-by-SNP cosine similarity matrix.
#' @param cluster Named integer vector of cluster assignments (`names` = SNP ids).
#' @param edge_threshold Absolute similarity used for connectedness and
#'   connected-component counts. Defaults to `0.5`.
#' @param min_module_size Minimum SNPs for `reliable` to be `TRUE`. Defaults to 3.
#' @param min_mean_internal Minimum mean internal similarity for `reliable`.
#'   Defaults to `0.3`.
#' @param min_connectedness Minimum pair-connectedness for `reliable`.
#'   Defaults to `0.5`.
#' @return A dataframe, one row per module, with internal/external similarity,
#'   connectedness, component coverage, diagnostic silhouette, and a `reliable`
#'   flag.
#' @export
summarise_snp_module_quality <- function(s_matrix,
                                         cluster,
                                         edge_threshold = 0.5,
                                         min_module_size = 3L,
                                         min_mean_internal = 0.3,
                                         min_connectedness = 0.5) {
  if (!is.matrix(s_matrix) || nrow(s_matrix) != ncol(s_matrix)) {
    stop("s_matrix must be a square matrix")
  }
  if (is.null(names(cluster))) {
    if (is.null(colnames(s_matrix))) {
      stop("cluster must be named, or s_matrix must have colnames")
    }
    names(cluster) <- colnames(s_matrix)
  }
  if (edge_threshold < 0) {
    stop("edge_threshold must be non-negative")
  }

  snp_ids <- intersect(names(cluster), colnames(s_matrix))
  if (length(snp_ids) < 2) {
    stop("Need at least 2 SNPs present in both cluster and s_matrix")
  }
  S <- s_matrix[snp_ids, snp_ids, drop = FALSE]
  diag(S) <- 1
  S[!is.finite(S)] <- 0
  cluster <- cluster[snp_ids]
  silhouette <- .silhouette_from_similarity(S, cluster)

  cluster_ids <- sort(unique(cluster))
  rows <- lapply(cluster_ids, function(cl) {
    members <- names(cluster)[cluster == cl]
    n <- length(members)
    rest <- setdiff(snp_ids, members)

    if (n >= 2) {
      internal <- S[members, members, drop = FALSE]
      internal_vals <- internal[upper.tri(internal)]
      mean_internal <- mean(internal_vals)
      n_edges <- sum(abs(internal_vals) >= edge_threshold)
      n_pairs <- length(internal_vals)
      connectedness <- n_edges / n_pairs
      binary <- abs(internal)
      diag(binary) <- 0
      binary[binary < edge_threshold] <- 0
      binary[binary > 0] <- 1
      graph <- igraph::graph_from_adjacency_matrix(
        binary,
        mode = "undirected",
        diag = FALSE
      )
      comp <- igraph::components(graph)
      n_components <- comp$no
      largest_frac <- max(comp$csize) / n
    } else {
      mean_internal <- NA_real_
      connectedness <- NA_real_
      n_edges <- 0L
      n_pairs <- 0L
      n_components <- 1L
      largest_frac <- 1
    }

    if (n >= 1 && length(rest) >= 1) {
      mean_external <- mean(S[members, rest, drop = FALSE])
    } else {
      mean_external <- NA_real_
    }

    mean_sil <- mean(silhouette[members], na.rm = TRUE)
    if (!is.finite(mean_sil)) {
      mean_sil <- NA_real_
    }

    reliable <- is.finite(n) && n >= min_module_size &&
      is.finite(mean_internal) && mean_internal >= min_mean_internal &&
      is.finite(connectedness) && connectedness >= min_connectedness

    return(data.frame(
      cluster = cl,
      n_snps = n,
      mean_internal_similarity = mean_internal,
      mean_external_similarity = mean_external,
      separation = mean_internal - mean_external,
      connectedness = connectedness,
      n_internal_edges = as.integer(n_edges),
      n_internal_pairs = as.integer(n_pairs),
      n_components = as.integer(n_components),
      largest_component_frac = largest_frac,
      mean_silhouette = mean_sil,
      reliable = reliable,
      stringsAsFactors = FALSE
    ))
  })

  result <- dplyr::bind_rows(rows) |>
    dplyr::arrange(dplyr::desc(reliable), dplyr::desc(n_snps))
  return(result)
}


.silhouette_from_similarity <- function(s_matrix, cluster) {
  snp_ids <- colnames(s_matrix)
  dist_mat <- 1 - s_matrix
  diag(dist_mat) <- 0
  sil <- stats::setNames(rep(NA_real_, length(snp_ids)), snp_ids)
  cluster_ids <- unique(cluster)

  for (id in snp_ids) {
    own <- as.character(cluster[[id]])
    members <- names(cluster)[as.character(cluster) == own]
    members <- setdiff(members, id)
    if (length(members) == 0) {
      next
    }
    a <- mean(dist_mat[id, members])
    other_means <- vapply(setdiff(as.character(cluster_ids), own), function(cl) {
      others <- names(cluster)[as.character(cluster) == cl]
      return(mean(dist_mat[id, others]))
    }, numeric(1))
    if (length(other_means) == 0 || !any(is.finite(other_means))) {
      next
    }
    b <- min(other_means, na.rm = TRUE)
    denom <- max(a, b)
    if (!is.finite(denom) || denom == 0) {
      sil[[id]] <- 0
    } else {
      sil[[id]] <- (b - a) / denom
    }
  }
  return(sil)
}
