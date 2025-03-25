library(ggplot2)
library(ggrepel)
library(dplyr)

plot_Manhattan = function(ordered_snp_data, sig_thres = c(0.05, 0.1, 0.2), plot_q = TRUE, path) {
  threshold_colors <- c("#4E79A7", "#E15759", "#76B7B2")
  thresholds <- data.frame(
    sig = sig_thres,
    color = threshold_colors[(seq_along(sig_thres)-1 %% length(threshold_colors)) + 1],
    linetype = "dashed"
  )
  
  # Data preprocessing ------------------------------------
  ordered_snp_data <- ordered_snp_data %>%
    mutate(chr = factor(chr, levels = unique(sort(as.numeric(chr)))))
  
  # Calculate Manhattan plot coordinates
  cumulative_pos <- function(chr, pos) {
    spacing <- 5e7  # Chromosome spacing: 50,000,000 bp
    df <- data.frame(chr = chr, pos = pos) %>%
      group_by(chr) %>%
      summarise(max_pos = max(pos)) %>%
      mutate(
        add_pos = cumsum(lag(max_pos + spacing, default = 0)),
        chr_center = add_pos + max_pos / 2
      )
    left_join(data.frame(chr = chr, pos = pos), df, by = "chr") %>%
      mutate(pos_cum = pos + add_pos) %>%
      pull(pos_cum)
  }
  
  # Add cumulative position and color grouping
  plot_data <- ordered_snp_data %>%
    arrange(as.numeric(chr), pos) %>%
    mutate(
      pos_cum = cumulative_pos(chr, pos),
      color_group = as.numeric(chr) %% 2,
      y_value = if (plot_q) -log10(q) else W
    )
  
  # Create chromosome label positions
  axis_df <- plot_data %>%
    group_by(chr) %>%
    summarize(center = mean(pos_cum))
  
  # Select significant SNPs for labeling
  label_data <- plot_data %>%
    filter(q <= (max(sig_thres)+0.01)) %>%
    group_by(gene.closer) %>%
    arrange(desc(W)) %>%
    slice_head(n = 1) %>%
    ungroup()

  
  # Create Manhattan plot ------------------------------------------------------------
  colors <- c("#8DA9DA", "#F2B2B3")
  
  manhattan_plot <- ggplot(plot_data, aes(x = pos_cum, y = y_value)) +
    geom_point(aes(color = factor(color_group)), size = 5, alpha = 0.7) +
    scale_color_manual(values = colors, guide = "none") +
    geom_text_repel(
      data = label_data,
      aes(label = gene.closer),
      size = 4.5,
      # fontface = "bold",
      min.segment.length = 0,
      point.padding = 0.2,
      box.padding = 0.5,
      max.overlaps = Inf,
      segment.color = "grey50"
    ) +
    scale_x_continuous(
      breaks = axis_df$center,
      labels = axis_df$chr,
      expand = expansion(mult = 0.02)
    ) +
    scale_y_continuous(expand = expansion(mult = c(0.05, 0.1))) +
    theme_classic(base_size = 12) +
    theme(
      legend.position = "none",
      panel.grid.major.y = element_line(color = "grey90", linewidth = 0.2),
      axis.line = element_line(color = "black", linewidth = 0.4),
      axis.ticks = element_line(color = "black", linewidth = 0.4),
      axis.text = element_text(color = "black", size = 14),
      axis.title = element_text(face = "bold", size = 16),
      plot.title = element_text(hjust = 0.5, face = "bold")
    ) +
    labs(
      x = "Chromosome",
      y = if (plot_q) expression(-log[10](q)) else "W"
    )
  
  # Add threshold lines and annotations
  if (plot_q) {
    for (i in 1:nrow(thresholds)) {
      manhattan_plot <- manhattan_plot +
        geom_hline(
          yintercept = -log10(thresholds$sig[i]),
          color = thresholds$color[i],
          linetype = thresholds$linetype[i],
          linewidth = 0.5
        ) +
        annotate("text",
                 x = -Inf, y = -log10(thresholds$sig[i]),
                 label = sprintf("FDR = %.2f", thresholds$sig[i]),
                 hjust = -0.1, vjust = -0.5,
                 size = 4, color = thresholds$color[i]
        )
    }
  } else {
    W_thresholds <- c()
    no_snps_indices <- c()
    for (i in 1:nrow(thresholds)) {
      sig_snps <- label_data$W[label_data$q <= thresholds$sig[i]+0.01]
      if (length(sig_snps) > 0) {
        W_thresholds <- c(W_thresholds, min(sig_snps))
      } else {
        W_thresholds <- c(W_thresholds, NA)
        no_snps_indices <- c(no_snps_indices, i)
      }
    }
    
    for (i in 1:length(W_thresholds)) {
      if (!is.na(W_thresholds[i])) {
        manhattan_plot <- manhattan_plot +
          geom_hline(
            yintercept = W_thresholds[i],
            color = thresholds$color[i],
            linetype = thresholds$linetype[i],
            linewidth = 0.5
          ) +
          annotate("text",
                   x = -Inf, y = W_thresholds[i],
                   label = sprintf("FDR = %.2f", thresholds$sig[i]),
                   hjust = -0.1, vjust = -0.5,
                   size = 4, color = thresholds$color[i]
          )
      }
    }
    
    if (length(no_snps_indices) > 0) {
      for (j in 1:length(no_snps_indices)) {
        i <- no_snps_indices[j]
        manhattan_plot <- manhattan_plot +
          annotate("text",
                   x = Inf, y = Inf,
                   label = sprintf("No SNPs reached FDR = %.2f", thresholds$sig[i]),
                   hjust = 1.1, vjust = 1.1 + (j-1)*0.5,
                   size = 4, color = thresholds$color[i]
          )
      }
    }
  }
  
  print(manhattan_plot)
  ggsave(path,
         plot = manhattan_plot,
         width = 12,
         height = 6,
         dpi = 900)
}



