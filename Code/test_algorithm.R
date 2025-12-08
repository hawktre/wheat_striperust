assign_predictions_to_sources <- function(error_matrix) {
  n_pred <- nrow(error_matrix)
  n_src <- ncol(error_matrix)
  
  # Initialize assignments
  assignments <- rep(NA, n_pred)
  names(assignments) <- rownames(error_matrix)
  assigned_sources <- logical(n_src)
  names(assigned_sources) <- colnames(error_matrix)
  
  # Step 1: Assign perfect matches (distance = 0)
  for (i in 1:n_pred) {
    for (j in 1:n_src) {
      if (error_matrix[i, j] == 0 && !assigned_sources[j]) {
        assignments[i] <- colnames(error_matrix)[j]
        assigned_sources[j] <- TRUE
        break
      }
    }
  }
  
  # Step 2: Assign remaining predictions iteratively
  while (any(is.na(assignments))) {
    unassigned <- which(is.na(assignments))
    available_sources <- which(!assigned_sources)
    
    if (length(available_sources) == 0) break
    
    # Find minimum distance among unassigned predictions to available sources
    min_dist <- Inf
    candidates <- list()
    
    for (i in unassigned) {
      for (j in available_sources) {
        dist <- error_matrix[i, j]
        if (dist < min_dist) {
          min_dist <- dist
          candidates <- list(list(pred = i, src = j))
        } else if (dist == min_dist) {
          candidates[[length(candidates) + 1]] <- list(pred = i, src = j)
        }
      }
    }
    
    # If only one candidate at minimum distance, assign it
    if (length(candidates) == 1) {
      pred_idx <- candidates[[1]]$pred
      src_idx <- candidates[[1]]$src
      assignments[pred_idx] <- colnames(error_matrix)[src_idx]
      assigned_sources[src_idx] <- TRUE
    } else {
      # Multiple candidates at minimum distance - need tie-breaking
      # Collect all unique predictions involved in ties
      pred_indices <- unique(sapply(candidates, function(x) x$pred))
      
      # For each tied prediction, find its next best alternative
      next_dists <- numeric(length(pred_indices))
      names(next_dists) <- pred_indices
      
      for (k in seq_along(pred_indices)) {
        pred_idx <- pred_indices[k]
        # Get distances to all available sources
        dists <- error_matrix[pred_idx, available_sources]
        # Remove the minimum distance (current tie)
        dists <- dists[dists > min_dist]
        next_dists[k] <- if (length(dists) > 0) min(dists) else Inf
      }
      
      # Sort predictions by worst alternatives (descending)
      sorted_order <- order(next_dists, decreasing = TRUE)
      sorted_preds <- pred_indices[sorted_order]
      
      # Assign only the first prediction (worst alternative) to its best tied source
      pred_to_assign <- sorted_preds[1]
      
      # Find which sources this prediction is tied for at min_dist
      tied_sources <- available_sources[error_matrix[pred_to_assign, available_sources] == min_dist]
      
      # Assign to the first tied source (they're all equally good at this distance)
      if (length(tied_sources) > 0) {
        assignments[pred_to_assign] <- colnames(error_matrix)[tied_sources[1]]
        assigned_sources[tied_sources[1]] <- TRUE
      }
      
      # Other predictions will be reconsidered in the next iteration
    }
  }
  
  # Create named vector with distances
  distances <- numeric(n_pred)
  names(distances) <- rownames(error_matrix)
  
  for (i in 1:n_pred) {
    if (!is.na(assignments[i])) {
      src_col <- which(colnames(error_matrix) == assignments[i])
      distances[i] <- error_matrix[i, src_col]
    } else {
      distances[i] <- NA
    }
  }
  
  return(distances)
}

# Example 1: Two-way tie (original example)
cat("=== Example 1: Two-way tie ===\n")
error <- matrix(c(
  0.0, 15.0, 20.0, 25.0,  # Prediction 1: clearly closest to source 1
  18.0,  0, 22.0, 30.0,   # Prediction 2: clearly closest to source 2
  10.0, 12.0, 10.0, 35.0, # Prediction 3: tied (10.0) between sources 1 and 3
  14.0, 16.0, 10.0, 40.0  # Prediction 4: also tied (10.0) for source 3
), nrow = 4, byrow = TRUE)

rownames(error) <- c("1", "2", "3", "4")
colnames(error) <- c("1", "2", "13", "9")

result <- assign_predictions_to_sources(error)
print("Distances to matched sources:")
print(result)

# Also show which source each prediction was matched to
cat("\nMatched sources:\n")
for (i in 1:length(result)) {
  # Find which source has this distance
  distances_for_pred <- error[i, ]
  matched_src <- names(distances_for_pred)[distances_for_pred == result[i]][1]
  cat("Prediction", names(result)[i], "-> Source", matched_src, 
      "(distance:", result[i], ")\n")
}

cat("\n")
cat("=== Example 2: Three-way tie ===\n")
error2 <- matrix(c(
  0.0, 20.0, 25.0, 30.0, 35.0,  # Prediction 1: perfect match to source 1
  15.0, 15.0, 15.0, 40.0, 45.0, # Prediction 2: 3-way tie at 15 (sources 1,2,3)
  15.0, 15.0, 15.0, 35.0, 50.0, # Prediction 3: 3-way tie at 15 (sources 1,2,3)
  15.0, 15.0, 15.0, 38.0, 42.0, # Prediction 4: 3-way tie at 15 (sources 1,2,3)
  22.0, 25.0, 28.0, 0.0, 55.0   # Prediction 5: perfect match to source 4
), nrow = 5, byrow = TRUE)

rownames(error2) <- c("A", "B", "C", "D", "E")
colnames(error2) <- c("1", "2", "3", "4", "5")

cat("\nPrediction 2, 3, 4 all tied at distance 15 for sources 1, 2, and 3\n")
cat("Next closest for Pred 2: source 4 at 40\n")
cat("Next closest for Pred 3: source 4 at 35\n")
cat("Next closest for Pred 4: source 4 at 38\n")
cat("\nExpected: Pred 2 should get one of the tied sources (worst alternative: 40)\n")
cat("Then Pred 4 (next worst: 38), then Pred 3 (best alternative: 35)\n\n")

result2 <- assign_predictions_to_sources(error2)
print("Distances to matched sources:")
print(result2)

cat("\nMatched sources:\n")
for (i in 1:length(result2)) {
  distances_for_pred <- error2[i, ]
  matched_src <- names(distances_for_pred)[distances_for_pred == result2[i]][1]
  cat("Prediction", names(result2)[i], "-> Source", matched_src, 
      "(distance:", result2[i], ")\n")
}