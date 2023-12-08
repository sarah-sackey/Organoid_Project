#geoms function
plot_col = function(df,x_name, y_name, col_by, geom_type, col_palette, fac_name, x_lab, y_lab, plot_title) {
  df |>
    ggplot(aes({{x_name}}, {{y_name}}, fill = col_by)) +
    geom_type + theme_bw() + theme(axis.text.x=element_text(angle=90, hjust=1)) +
    col_palette + facet_grid(as.formula(paste("~", fac_name))) + xlab(x_lab) + ylab(y_lab) + ggtitle(plot_title)
}
plot_col(df = CD4_datalong, x_name = factor(p24_type, levels = c("Uninfected_CD4Total", "Infected_p24neg", "Infected_p24pos")), 
         y_name = freq, col_by = CD4_datalong$cell_type2, geom_type = geom_bar(position = "fill", stat= "identity"), 
         col_palette = scale_fill_moma_d("Smith"), fac_name = "Timepoint", x_lab = "MOI", y_lab = "% of CD4 subset",
         plot_title = "CD4 cell subset over time")
