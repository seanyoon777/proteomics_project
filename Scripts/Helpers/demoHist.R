
generate_sampleHist <- function(patientdata, drawdate_windows, ymin = 30, ymax = 100, yint = 10) {
 
  drawdate_hists <- list()
  for (i in 1:(length(drawdate_windows) - 1)) {
    patientdata_bin <- patientdata[
      abs(patientdata$drawdate_diff) < drawdate_windows[i + 1] & 
        abs(patientdata$drawdate_diff) >= drawdate_windows[i], ]
    
    drawdate_hists[[i]] <- patientdata_bin %>% 
      ggplot(aes(x = avg_drawage, fill = final_status)) + 
      geom_histogram(binwidth = 1, breaks = seq(ymin, ymax, by = yint), 
                     position = "stack", show.legend = FALSE, color = "black") + 
      labs(title = paste0("Drawn within ", drawdate_windows[i], " and ", 
                          drawdate_windows[i + 1], " days"), 
           x = "Average age when drawn", y = "Frequency") +
      theme_classic() +
      theme(panel.border = element_rect(colour = "black", fill = NA, size = 1),
            plot.background = element_rect(fill = "white"),
            axis.title = element_text(size = 9.5),
            axis.text = element_text(size = 7),
            plot.title = element_text(size = 10, hjust = 0.5), 
            axis.line = element_blank()) +
      scale_y_continuous(breaks = seq(0, 100, 10)) +
      scale_fill_manual(values = c("CO" = "#00BFC4", "AD" = "#F8766D"))
    
    maxcount <- max(hist(patientdata_bin$avg_drawage, 
                         breaks = seq(30, 100, by = 10))$counts)
    drawdate_hists[[i]] <- drawdate_hists[[i]] + 
      annotate("text", x = 95, y = maxcount, label = paste0("n = ", nrow(patientdata_bin)), 
               size = 3)
  }
  #AD = red, CO = blue
  gridExtra::grid.arrange(grobs = drawdate_hists, ncol = 3)
}
