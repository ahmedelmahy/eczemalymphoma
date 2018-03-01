require(FirebrowseR)
# get projects names
Metadata.Cohorts(format = "csv")

# iterate over the pages recieved from the api to get the whole dataset
# for project skin cutaneous melanoma, we use:
all.Received = F
page.Counter = 1
page.size = 150
skcm.Pats = list()
while(all.Received == F){
    skcm.Pats[[page.Counter]] = Samples.Clinical(format = "csv",
                                                 cohort = "SKCM",
                                                 page_size = page.size,
                                                 page = page.Counter)
    if(page.Counter > 1)
        colnames(skcm.Pats[[page.Counter]]) = colnames(skcm.Pats[[page.Counter-1]])

    if(nrow(skcm.Pats[[page.Counter]]) < page.size){
        all.Received = T
    } else{
        page.Counter = page.Counter + 1
    }
}
skcm.Pats = do.call(rbind, skcm.Pats)
dim(skcm.Pats)
colnames(skcm.Pats)
skcm.Pats$person_neoplasm_cancer_status
