source("src/PDPD.R")
kernel_dir ="src/Data/Kernels"
kernel_Rdata_dir ="src/Data/Kernels/RData"
kernel_dir_already_exists = dir.exists(kernel_dir)
cat("\nKernel folder already exists :",kernel_dir_already_exists)
if(!kernel_dir_already_exists){
    cat("\nKernel folder creating....")
    dir.create(kernel_dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
    dir.create(kernel_Rdata_dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
}
main()
