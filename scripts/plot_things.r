args <- commandArgs(trailingOnly = TRUE)

setwd(getwd())

data <- read.csv(args[1], header = TRUE, sep = ',')

#write.csv(data, file = args[2], row.names = FALSE)
pdf(file = args[2])
plot(data$bill_length_mm, data$bill_depth_mm, pch = 19, col = "blue",
     xlab = "Bill Length (mm)", ylab = "Bill Depth (mm)",
     main = "Penguins: Bill Length vs. Bill Depth")
dev.off()

