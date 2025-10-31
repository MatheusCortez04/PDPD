source("src/PDPD.R")
kernel_dir ="src/Data/Kernels"
kernel_Rdata_dir ="src/Data/Kernels/RData"
kernel_dir_already_exists = dir.exists(kernel_dir)

drug_dir ="src/Data/Drug"
drug_Rdata_dir ="src/Data/Drug/RData"
drug_dir_already_exists = dir.exists(drug_dir)
cat("\nKernel folder already exists :",kernel_dir_already_exists)
cat("\nDrug folder already exists :",drug_dir_already_exists)
if(!kernel_dir_already_exists){
    cat("\nKernel folder creating....")
    dir.create(kernel_dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
    dir.create(kernel_Rdata_dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
}
if(!drug_dir_already_exists){
    cat("\nKernel folder creating....")
    dir.create(drug_dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
    dir.create(drug_Rdata_dir, showWarnings = TRUE, recursive = FALSE, mode = "0777")
}
main()
