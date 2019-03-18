#------------------------------------------#
# DNA methylation: Get all TCGA IDAT files #
#------------------------------------------#
projects <- TCGAbiolinks:::getGDCprojects()$project_id
projects <- projects[grepl('^TCGA',projects,perl=T)]
match.file.cases.all <- NULL
for(proj in projects){
print(proj)
query <- GDCquery(project = proj,
data.category = "Raw microarray data",
data.type = "Raw intensities",
experimental.strategy = "Methylation array",
legacy = TRUE,
file.type = ".idat",
platform = "Illumina Human Methylation 450")
match.file.cases <- getResults(query,cols=c("cases","file_name"))
match.file.cases$project <- proj
match.file.cases.all <- rbind(match.file.cases.all,match.file.cases)
tryCatch(GDCdownload(query, method = "api", files.per.chunk = 20),
error = function(e) GDCdownload(query, method = "client"))
}
# This will create a map between idat file name, cases (barcode) and project
readr::write_tsv(match.file.cases.all, path =  "idat_filename_case.txt")
# code to move all files to local folder
for(file in dir(".",pattern = ".idat", recursive = T)){
TCGAbiolinks::move(file,basename(file))
}


#------------------------------------------#
# DNA methylation: aligned against hg19    #
#------------------------------------------#
query_meth.hg19 <- GDCquery(project= "TCGA-READ",
data.category = "DNA methylation",
platform = "Illumina Human Methylation 450",
barcode = c("TCGA-DY-A1DD-01A-21D-A153-05",
"TCGA-AG-A020-11A-11D-A081-05","TCGA-AG-3591-01A-01D-1734-05","TCGA-AG-A01W-01A-21D-A081-05","TCGA-AG-A02N-11A-11D-A081-05","TCGA-EF-5831-01A-01D-1658-05",
"TCGA-EI-6507-01A-11D-1734-05","TCGA-AG-3725-01A-11D-1734-05","TCGA-AH-6549-01A-11D-1828-05","TCGA-G5-6572-02A-12D-1828-05","TCGA-CI-6622-01A-11D-1828-05",
"TCGA-AG-A02N-01A-11D-A081-05","TCGA-DC-6683-01A-11D-1828-05","TCGA-EI-6511-01A-11D-1734-05","TCGA-BM-6198-01A-11D-1734-05","TCGA-EI-6917-01A-11D-1926-05",
"TCGA-AG-3732-01A-11D-1658-05","TCGA-F5-6702-01A-11D-1828-05","TCGA-AG-3725-11A-01D-1734-05","TCGA-EI-6882-01A-11D-1926-05","TCGA-DY-A1DC-01A-31D-A153-05",
"TCGA-DY-A0XA-01A-11D-A153-05","TCGA-AF-A56K-01A-32D-A39G-05","TCGA-CL-5918-01A-11D-1658-05","TCGA-EI-6884-01A-11D-1926-05","TCGA-EI-6885-01A-11D-1926-05",
"TCGA-EI-6509-01A-11D-1734-05","TCGA-AF-A56N-01A-12D-A39G-05","TCGA-CI-6621-01A-11D-1828-05","TCGA-AG-A01Y-11A-11D-A081-05","TCGA-AG-A01Y-01A-41D-A081-05",
"TCGA-AH-6903-01A-11D-1926-05","TCGA-AF-3911-01A-01D-1734-05","TCGA-DC-5337-01A-01D-1658-05","TCGA-DY-A1DG-01A-11D-A153-05","TCGA-AG-A036-01A-12D-A081-05",
"TCGA-G5-6641-01A-11D-1828-05","TCGA-EI-7004-01A-11D-1926-05","TCGA-DC-6158-01A-11D-1658-05","TCGA-AG-A01W-11A-11D-A081-05","TCGA-AF-2693-01A-02D-1734-05",
"TCGA-DC-4745-01A-01D-1734-05","TCGA-F5-6464-01A-11D-1734-05","TCGA-DC-4749-01A-01D-1734-05","TCGA-AG-A020-01A-21D-A081-05","TCGA-AH-6643-01A-11D-1828-05",
"TCGA-EI-6514-01A-11D-1734-05","TCGA-AH-6544-01A-11D-1828-05","TCGA-EF-5830-01A-01D-1658-05","TCGA-DT-5265-01A-21D-1828-05","TCGA-AG-4022-01A-01D-1734-05",
"TCGA-CL-5917-01A-11D-1658-05","TCGA-EI-6513-01A-21D-1734-05","TCGA-F5-6861-01A-11D-1926-05","TCGA-G5-6572-01A-11D-1828-05","TCGA-CI-6624-01C-11D-1828-05",
"TCGA-AG-A036-11A-11D-A081-05","TCGA-DC-6160-01A-11D-1658-05","TCGA-CL-4957-01A-01D-1734-05","TCGA-AG-A026-01A-32D-A081-05","TCGA-F5-6863-01A-11D-1926-05",
"TCGA-DC-6154-01A-31D-1926-05","TCGA-F5-6814-01A-31D-1926-05","TCGA-DY-A1DF-01A-11D-A153-05","TCGA-AG-4021-01A-01D-1734-05","TCGA-F5-6864-01A-11D-1926-05",
"TCGA-DC-6682-01A-11D-1828-05","TCGA-EI-6881-01A-11D-1926-05","TCGA-AH-6547-01A-11D-1828-05","TCGA-AH-6897-01A-11D-1926-05","TCGA-CI-6623-01B-11D-1828-05",
"TCGA-DC-5869-01A-01D-1658-05","TCGA-F5-6571-01A-12D-1828-05","TCGA-AG-3731-01A-11D-1734-05","TCGA-EI-6883-01A-31D-1926-05","TCGA-AH-6644-01A-21D-1828-05",
"TCGA-EI-6506-01A-11D-1734-05","TCGA-G5-6233-01A-11D-1734-05","TCGA-DY-A1DE-01A-11D-A153-05","TCGA-DC-6157-01A-11D-1658-05","TCGA-AF-6136-01A-11D-1828-05",
"TCGA-EI-6512-01A-11D-1734-05","TCGA-F5-6810-01A-11D-1828-05","TCGA-EI-6510-01A-11D-1734-05","TCGA-DC-6681-01A-11D-1828-05","TCGA-F5-6465-01A-11D-1734-05",
"TCGA-AG-3731-11A-01D-1734-05","TCGA-AF-2687-01A-02D-1734-05","TCGA-DC-6156-01A-11D-1658-05","TCGA-AG-3742-01A-11D-1658-05","TCGA-AF-2690-01A-02D-1734-05",
"TCGA-EI-6508-01A-11D-1734-05","TCGA-AF-6655-01A-11D-1828-05","TCGA-AG-3592-01A-02D-1734-05","TCGA-AF-6672-01A-11D-1828-05","TCGA-CI-6620-01A-11D-1828-05",
"TCGA-CI-6619-01B-11D-1828-05","TCGA-DY-A1H8-01A-21D-A153-05","TCGA-AF-A56L-01A-31D-A39G-05","TCGA-F5-6812-01A-11D-1828-05","TCGA-EI-7002-01A-11D-1926-05",
"TCGA-DC-6155-01A-11D-1658-05","TCGA-AF-4110-01A-02D-1734-05","TCGA-G5-6235-01A-11D-1734-05","TCGA-F5-6813-01A-11D-1828-05","TCGA-F5-6811-01A-11D-1828-05"),
legacy = TRUE)
GDCdownload(query_meth.hg19)
data.hg19 <- GDCprepare(query_meth.hg19)









"TCGA-AA-3712-01A-21D-1721-05","TCGA-CK-6747-01A-11D-1837-05","TCGA-AA-3502-11A-01D-1407-05","TCGA-D5-6536-01A-11D-1721-05","TCGA-CM-6676-01A-11D-1837-05",
"TCGA-G4-6306-01A-11D-1772-05","TCGA-AZ-6605-01A-11D-1837-05","TCGA-G4-6588-01A-11D-1772-05","TCGA-G4-6297-01A-11D-1721-05","TCGA-A6-2675-01A-02D-1721-05",
"TCGA-DM-A0X9-01A-11D-A153-05","TCGA-AA-3697-01A-01D-1721-05","TCGA-G4-6315-01A-11D-1721-05","TCGA-F4-6808-01A-11D-1837-05","TCGA-CM-5864-01A-01D-1651-05",
"TCGA-G4-6322-01A-11D-1721-05","TCGA-DM-A1DB-01A-11D-A153-05","TCGA-NH-A6GB-01A-11D-A36Y-05","TCGA-AA-3662-01A-01D-1721-05","TCGA-AZ-6608-01A-11D-1837-05",
"TCGA-A6-5661-01A-01D-1651-05","TCGA-A6-6780-01A-11D-A27A-05","TCGA-A6-5660-01A-01D-1651-05","TCGA-AZ-4308-01A-01D-1407-05","TCGA-D5-6534-01A-21D-1926-05",
"TCGA-G4-6627-01A-11D-1772-05","TCGA-CM-6677-01A-11D-1837-05","TCGA-SS-A7HO-01A-21D-A36Y-05","TCGA-D5-5541-01A-01D-1651-05","TCGA-F4-6806-01A-11D-1837-05",
"TCGA-A6-2686-11A-01D-1551-05","TCGA-A6-2671-01A-01D-1407-05","TCGA-D5-6930-01A-11D-1926-05","TCGA-A6-A566-01A-11D-A28O-05","TCGA-A6-5662-01A-01D-1651-05",
"TCGA-G4-6309-01A-21D-1837-05","TCGA-F4-6463-01A-11D-1721-05","TCGA-DM-A28M-01A-12D-A16X-05","TCGA-D5-6927-01A-21D-1926-05","TCGA-AA-3492-01A-01D-1407-05",
"TCGA-CK-4951-01A-01D-1407-05","TCGA-A6-2685-01A-01D-1407-05","TCGA-A6-2684-01A-01D-1407-05","TCGA-A6-5667-11A-01D-1721-05","TCGA-AA-3663-11A-01D-1721-05",
"TCGA-D5-6924-01A-11D-1926-05","TCGA-AA-3510-11A-01D-1407-05","TCGA-A6-4107-01A-02D-1407-05","TCGA-A6-5667-01A-21D-1721-05","TCGA-AY-A71X-01A-12D-A36Y-05",
"TCGA-DM-A28A-01A-21D-A16X-05","TCGA-A6-6649-01A-11D-1772-05","TCGA-D5-5538-01A-01D-1651-05","TCGA-D5-6926-01A-11D-1926-05","TCGA-G4-6321-01A-11D-1721-05",
"TCGA-AZ-4615-01A-01D-1407-05","TCGA-A6-6648-01A-11D-1772-05","TCGA-A6-6652-01A-11D-1772-05","TCGA-F4-6855-01A-11D-1926-05","TCGA-AZ-6601-01A-11D-1772-05",
"TCGA-D5-6531-01A-11D-1721-05","TCGA-AD-6888-01A-11D-1926-05","TCGA-D5-6920-01A-11D-1926-05","TCGA-CA-6716-01A-11D-1837-05","TCGA-AZ-4614-01A-01D-1407-05",
"TCGA-5M-AATE-01A-11D-A40X-05","TCGA-CK-6748-01A-11D-1837-05","TCGA-D5-6537-01A-11D-1721-05","TCGA-AA-3506-11A-01D-1407-05","TCGA-3L-AA1B-01A-11D-A36Y-05",
"TCGA-AA-3495-11A-01D-1407-05","TCGA-G4-6320-01A-11D-1721-05","TCGA-AZ-6603-01A-11D-1837-05","TCGA-AZ-6601-11A-01D-1772-05","TCGA-D5-6538-01A-11D-1721-05",
"TCGA-DM-A282-01A-12D-A16X-05","TCGA-A6-5659-01A-01D-1651-05","TCGA-CK-6746-01A-11D-1837-05","TCGA-AA-3509-11A-01D-1407-05","TCGA-A6-2677-01A-01D-A27A-05",
"TCGA-A6-6780-01A-11D-1837-05","TCGA-F4-6854-01A-11D-1926-05","TCGA-AM-5821-01A-01D-1651-05","TCGA-AZ-5403-01A-01D-1651-05","TCGA-A6-6781-01A-22D-A27A-05",
"TCGA-A6-5666-01A-01D-1651-05","TCGA-DM-A1D9-01A-11D-A153-05","TCGA-CA-6718-01A-11D-1837-05","TCGA-QG-A5Z2-01A-11D-A28O-05","TCGA-G4-6311-01A-11D-1721-05",
"TCGA-A6-6650-01A-11D-1772-05","TCGA-A6-5661-01B-05D-2299-05","TCGA-AA-3494-11A-01D-1407-05","TCGA-CM-6171-01A-11D-1651-05","TCGA-D5-6541-01A-11D-1721-05",
"TCGA-CA-6715-01A-21D-1837-05","TCGA-CK-4948-01B-11D-1651-05","TCGA-G4-6298-11A-01D-1721-05","TCGA-CA-6717-01A-11D-1837-05","TCGA-CM-6162-01A-11D-1651-05",
"TCGA-CM-5862-01A-01D-1651-05","TCGA-D5-6923-01A-11D-1926-05","TCGA-AZ-6599-11A-01D-1772-05","TCGA-D5-6931-01A-11D-1926-05","TCGA-F4-6703-01A-11D-1837-05",
"TCGA-DM-A1D0-01A-11D-A153-05","TCGA-G4-6314-01A-11D-1721-05","TCGA-AU-3779-01A-01D-1721-05","TCGA-D5-5540-01A-01D-1651-05","TCGA-DM-A285-01A-11D-A16X-05",
"TCGA-AA-3496-01A-21D-1837-05","TCGA-CK-5914-01A-11D-1651-05","TCGA-WS-AB45-01A-11D-A40X-05","TCGA-DM-A1DA-01A-11D-A153-05","TCGA-AY-A8YK-01A-11D-A40X-05",
"TCGA-AM-5820-01A-01D-1651-05","TCGA-A6-4105-01A-02D-1772-05","TCGA-NH-A50V-01A-11D-A28O-05","TCGA-A6-2682-11A-01D-1551-05","TCGA-AA-3660-01A-01D-1721-05",
"TCGA-DM-A1D8-01A-11D-A153-05","TCGA-AA-3494-01A-01D-1407-05","TCGA-AZ-4681-01A-01D-1407-05","TCGA-A6-6781-01B-06D-A27A-05","TCGA-AA-3502-01A-01D-1407-05",
"TCGA-D5-7000-01A-11D-1926-05","TCGA-A6-4107-11A-01D-1551-05","TCGA-F4-6460-01A-11D-1772-05","TCGA-DM-A28K-01A-21D-A16X-05","TCGA-G4-6293-01A-11D-1721-05",
"TCGA-A6-6782-01A-11D-1837-05","TCGA-CM-6168-01A-11D-1651-05","TCGA-A6-2680-01A-01D-1407-05","TCGA-NH-A50U-01A-33D-A36Y-05","TCGA-A6-A5ZU-01A-11D-A28O-05",
"TCGA-CK-5915-01A-11D-1651-05","TCGA-AD-6963-01A-11D-1926-05","TCGA-NH-A6GC-01A-12D-A40X-05","TCGA-AA-3655-01A-02D-1721-05","TCGA-CM-4744-01A-01D-1407-05",
"TCGA-CM-4748-01A-01D-1407-05","TCGA-A6-5657-01A-01D-1651-05","TCGA-D5-5539-01A-01D-1651-05","TCGA-A6-2679-01A-02D-1407-05","TCGA-CK-5913-01A-11D-1651-05",
"TCGA-DM-A1HB-01A-22D-A17Z-05","TCGA-A6-A56B-01A-31D-A28O-05","TCGA-DM-A1D6-01A-21D-A153-05","TCGA-A6-6141-01A-11D-1772-05","TCGA-AZ-6598-01A-11D-1772-05",
"TCGA-G4-6626-01A-11D-1772-05","TCGA-G4-6317-02A-11D-2064-05","TCGA-AZ-6599-01A-11D-1772-05","TCGA-AZ-6598-11A-01D-1772-05","TCGA-A6-5656-01A-21D-A27A-05",
"TCGA-D5-6539-01A-11D-1721-05","TCGA-CM-6161-01A-11D-1651-05","TCGA-D5-6532-01A-11D-1721-05","TCGA-G4-6298-01A-11D-1721-05","TCGA-CA-6719-01A-11D-1837-05",
"TCGA-CM-4746-01A-01D-1407-05","TCGA-DM-A1D7-01A-11D-A153-05","TCGA-NH-A8F7-01A-11D-A40X-05","TCGA-CM-6163-01A-11D-1651-05","TCGA-DM-A1HA-01A-11D-A153-05",
"TCGA-A6-2681-01A-01D-1407-05","TCGA-F4-6461-01A-11D-1772-05","TCGA-G4-6314-11A-01D-1721-05","TCGA-A6-2671-11A-01D-1551-05","TCGA-A6-A565-01A-31D-A28O-05",
"TCGA-AZ-6606-01A-11D-1837-05","TCGA-A6-2679-11A-01D-1551-05","TCGA-CM-5344-01A-21D-1721-05","TCGA-AA-3697-11A-01D-1721-05","TCGA-G4-6307-01A-11D-1721-05",
"TCGA-AA-3488-11A-01D-1407-05","TCGA-G4-6297-11A-01D-1721-05","TCGA-G4-6320-11A-01D-1721-05","TCGA-A6-6653-01A-11D-1772-05","TCGA-NH-A8F8-01A-72D-A40X-05",
"TCGA-4N-A93T-01A-11D-A36Y-05","TCGA-QG-A5YV-01A-11D-A28O-05","TCGA-DM-A28G-01A-11D-A16X-05","TCGA-AY-A54L-01A-11D-A28O-05","TCGA-AA-3495-01A-01D-1407-05",
"TCGA-AY-6386-01A-21D-1721-05","TCGA-CK-5912-01A-11D-1651-05","TCGA-NH-A8F7-06A-31D-A40X-05","TCGA-AZ-5407-01A-01D-1721-05","TCGA-A6-3810-01A-01D-A27A-05",
"TCGA-A6-2672-01B-03D-2299-05","TCGA-AZ-6607-01A-11D-1837-05","TCGA-AD-A5EJ-01A-11D-A28O-05","TCGA-5M-AAT4-01A-11D-A40X-05","TCGA-A6-5665-01B-03D-2299-05",
"TCGA-AA-3713-11A-01D-1721-05","TCGA-G4-6299-01A-11D-1772-05","TCGA-CA-5255-01A-11D-1837-05","TCGA-CM-5861-01A-01D-1651-05","TCGA-QL-A97D-01A-12D-A40X-05",
"TCGA-A6-2677-01B-02D-A27A-05","TCGA-DM-A288-01A-11D-A16X-05","TCGA-F4-6569-01A-11D-1772-05","TCGA-A6-6650-01A-11D-A27A-05","TCGA-CM-6169-01A-11D-1651-05",
"TCGA-AY-5543-01A-01D-1651-05","TCGA-RU-A8FL-01A-11D-A36Y-05","TCGA-CM-4750-01A-01D-1407-05","TCGA-G4-6322-11A-01D-1721-05","TCGA-A6-A567-01A-31D-A28O-05",
"TCGA-CM-4747-01A-01D-1407-05","TCGA-A6-5664-01A-21D-1837-05","TCGA-NH-A50T-01A-11D-A28O-05","TCGA-CK-4947-01B-11D-1651-05","TCGA-G4-6294-01A-11D-1772-05",
"TCGA-CM-5868-01A-01D-1651-05","TCGA-A6-6651-01A-21D-1837-05","TCGA-AA-3506-01A-01D-1407-05","TCGA-G4-6323-01A-11D-1721-05","TCGA-T9-A92H-01A-11D-A36Y-05",
"TCGA-AZ-4684-01A-01D-1407-05","TCGA-AD-6965-01A-11D-1926-05","TCGA-G4-6295-11A-01D-1721-05","TCGA-QG-A5YW-01A-11D-A28O-05","TCGA-G4-6295-01A-11D-1721-05",
"TCGA-G4-6628-01A-11D-1837-05","TCGA-CM-5341-01A-01D-1407-05","TCGA-CM-6166-01A-11D-1651-05","TCGA-A6-2675-11A-01D-1721-05","TCGA-D5-6529-01A-11D-1772-05",
"TCGA-CA-5256-01A-01D-1407-05","TCGA-AD-6889-01A-11D-1926-05","TCGA-D5-6929-01A-31D-1926-05","TCGA-AA-3663-01A-01D-1721-05","TCGA-A6-2684-01C-08D-A27A-05",
"TCGA-A6-5656-01B-02D-A27A-05","TCGA-D5-6530-01A-11D-1721-05","TCGA-G4-6586-01A-11D-1772-05","TCGA-F4-6856-01A-11D-1926-05","TCGA-CM-4752-01A-01D-1407-05",
"TCGA-CM-6165-01A-11D-1651-05","TCGA-DM-A0XD-01A-12D-A153-05","TCGA-AD-6901-01A-11D-1926-05","TCGA-CK-6751-01A-11D-1837-05","TCGA-AZ-6600-01A-11D-1772-05",
"TCGA-CM-4743-01A-01D-1721-05","TCGA-AD-6964-01A-11D-1926-05","TCGA-NH-A5IV-01A-42D-A36Y-05","TCGA-AD-A5EK-01A-11D-A28O-05","TCGA-G4-6625-11A-01D-1772-05",
"TCGA-F4-6704-01A-11D-1837-05","TCGA-AA-3713-01A-21D-1721-05","TCGA-F4-6805-01A-11D-1837-05","TCGA-5M-AAT5-01A-21D-A40X-05","TCGA-CA-5797-01A-01D-1651-05",
"TCGA-A6-3809-01B-04D-A27A-05","TCGA-AD-6890-01A-11D-1926-05","TCGA-F4-6809-01A-11D-1837-05","TCGA-QG-A5YX-01A-11D-A28O-05","TCGA-D5-6922-01A-11D-1926-05",
"TCGA-A6-2684-11A-01D-1551-05","TCGA-A6-2682-01A-01D-1407-05","TCGA-AZ-4616-01A-21D-1837-05","TCGA-AY-6197-01A-11D-1721-05","TCGA-G4-6625-01A-21D-1772-05",
"TCGA-CM-6680-01A-11D-1837-05","TCGA-CM-5348-01A-21D-1721-05","TCGA-A6-5665-01A-01D-1651-05","TCGA-A6-2681-11A-01D-1551-05","TCGA-AZ-4313-01A-01D-1407-05",
"TCGA-CM-6678-01A-11D-1837-05","TCGA-AA-3510-01A-01D-1407-05","TCGA-G4-6317-01A-11D-1721-05","TCGA-CM-6675-01A-11D-1837-05","TCGA-A6-6142-01A-11D-1772-05",
"TCGA-A6-6650-01B-02D-A27A-05","TCGA-CM-6674-01A-11D-1837-05","TCGA-CK-4950-01A-01D-1721-05","TCGA-D5-6928-01A-11D-1926-05","TCGA-A6-6781-01A-22D-1926-05",
"TCGA-A6-3810-01B-04D-A27A-05","TCGA-AD-6899-01A-11D-1926-05","TCGA-CM-5349-01A-21D-1721-05","TCGA-A6-5656-01A-21D-1837-05","TCGA-A6-2686-01A-01D-1407-05",
"TCGA-F4-6570-01A-11D-1772-05","TCGA-CM-6172-01A-11D-1651-05","TCGA-A6-6780-01B-04D-A27A-05","TCGA-CM-4751-01A-02D-1837-05","TCGA-A6-2685-11A-01D-1551-05",
"TCGA-G4-6310-01A-11D-1721-05","TCGA-AD-6548-01A-11D-1837-05","TCGA-F4-6807-01A-11D-1837-05","TCGA-CA-5796-01A-01D-1651-05","TCGA-CM-6679-01A-11D-1837-05",
"TCGA-CK-4952-01A-01D-1721-05","TCGA-G4-6302-11A-01D-1721-05","TCGA-CM-5860-01A-01D-1651-05","TCGA-F4-6459-01A-11D-1772-05","TCGA-DM-A28H-01A-11D-A16X-05",
"TCGA-AA-3492-11A-01D-1407-05","TCGA-4T-AA8H-01A-11D-A40X-05","TCGA-AY-A69D-01A-11D-A36Y-05","TCGA-A6-6138-01A-11D-1772-05","TCGA-A6-2680-11A-01D-1551-05",
"TCGA-D5-6535-01A-11D-1721-05","TCGA-AA-3511-01A-21D-1837-05","TCGA-CM-6164-01A-11D-1651-05","TCGA-AZ-4682-01B-01D-1407-05","TCGA-DM-A28C-01A-11D-A16X-05",
"TCGA-AA-3712-11A-01D-1721-05","TCGA-G4-6304-01A-11D-1926-05","TCGA-A6-6140-01A-11D-1772-05","TCGA-CM-6167-01A-11D-1651-05","TCGA-5M-AATA-01A-31D-A40X-05",
"TCGA-DM-A0XF-01A-11D-A153-05","TCGA-AZ-6600-11A-01D-1772-05","TCGA-AA-3489-01A-21D-1837-05","TCGA-CK-5916-01A-11D-1651-05","TCGA-AY-6196-01A-11D-1721-05",
"TCGA-A6-6137-01A-11D-1772-05","TCGA-AZ-4323-01A-21D-1837-05","TCGA-AA-3655-11A-01D-1721-05","TCGA-NH-A6GA-01A-11D-A36Y-05","TCGA-G4-6303-01A-11D-1772-05",
"TCGA-AU-6004-01A-11D-1721-05","TCGA-D5-5537-01A-21D-1926-05","TCGA-D5-6533-01A-11D-1721-05","TCGA-DM-A280-01A-12D-A16X-05","TCGA-AA-3488-01A-01D-1407-05",
"TCGA-AZ-4315-01A-01D-1407-05","TCGA-CA-5254-01A-21D-1837-05","TCGA-A6-6654-01A-21D-1837-05","TCGA-CM-6170-01A-11D-1651-05","TCGA-QG-A5Z1-01A-11D-A28O-05",
"TCGA-DM-A28F-01A-11D-A16X-05","TCGA-G4-6302-01A-11D-1721-05","TCGA-DM-A1D4-01A-21D-A153-05","TCGA-D5-6932-01A-11D-1926-05","TCGA-A6-5659-01A-01D-A27A-05",
"TCGA-AD-5900-01A-11D-1651-05","TCGA-5M-AAT6-01A-11D-A40X-05","TCGA-A6-5659-01B-04D-A27A-05","TCGA-DM-A28E-01A-11D-A16X-05","TCGA-AD-6895-01A-11D-1926-05",
"TCGA-A6-3809-01A-01D-A27A-05")




#------------------------#
#   Protein expression   #
#------------------------#
query <- GDCquery(project = "TCGA-READ",
data.category = "Protein expression",
legacy = TRUE,
barcode = c("TCGA-OX-A56R-01A-21-A44T-20","TCGA-08-0357-01A-21-1898-20"))
GDCdownload(query)
data <- GDCprepare(query, save = TRUE,
save.filename = "readProteinExpression.rda",
remove.files.prepared = TRUE)


#-----------------------------------------#
# Gene expression: aligned against hg19   #
#-----------------------------------------#
query.exp.hg19 <- GDCquery(project = "TCGA-READ",
data.category = "Gene expression",
data.type = "Gene expression quantification",
platform = "Illumina HiSeq",
file.type  = "normalized_results",
experimental.strategy = "RNA-Seq",
barcode = c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"),
legacy = TRUE)
GDCdownload(query.exp.hg19)
data <- GDCprepare(query.exp.hg19)


#-----------------------------------------#
#    Harmonized database | Copy Number    #
#-----------------------------------------#
query_READ_copy_segment <- GDCquery(project = "TCGA-READ",
data.category = "Copy Number Variation",
data.type = "Copy Number Segment",
barcode = c( "TCGA-OR-A5KU-01A-11D-A29H-01", "TCGA-OR-A5JK-01A-11D-A29H-01"))
GDCdownload(query_READ_copy_segment)
data_READ_copy_segment <- GDCprepare(query_READ_copy_segment)


query_COAD_copy_segment <- GDCquery(project = "TCGA-COAD",
data.category = "Copy Number Variation",
data.type = "Copy Number Segment",
barcode = c( "TCGA-OR-A5KU-01A-11D-A29H-01", "TCGA-OR-A5JK-01A-11D-A29H-01"))
GDCdownload(query_COAD_copy_segment)
data <- GDCprepare(query_READ_copy_segment)




#----------------#
#    GISTIC2     #
#----------------#
query_READ_GISTIC2 <- GDCquery(project = "TCGA-READ",
data.category = "Copy Number Variation",
data.type = "Gene Level Copy Number Scores",
access="open")
GDCdownload(query_READ_GISTIC2)
data_READ_GISTIC2 <- GDCprepare(query_READ_GISTIC2)



data_COAD_GISTIC2 <- GDCquery(project = "TCGA-COAD",
data.category = "Copy Number Variation",
data.type = "Gene Level Copy Number Scores",
access="open")
GDCdownload(data_COAD_GISTIC2)
data_COAD_GISTIC2 <- GDCprepare(query_COAD_GISTIC2)



#--------------------------------------------#
#   Gene expression: aligned against hg38    #
#--------------------------------------------#

# mRNA pipeline: https://gdc-docs.nci.nih.gov/Data/Bioinformatics_Pipelines/Expression_mRNA_Pipeline/
query.READ_exp.hg38 <- GDCquery(project = "TCGA-READ",
data.category = "Transcriptome Profiling",
data.type = "Gene Expression Quantification",
workflow.type = "HTSeq - FPKM-UQ",
barcode =  c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"))
GDCdownload(query.READ_exp.hg38)
expdat <- GDCprepare(query = query.READ_exp.hg38,
save = TRUE,
save.filename = "READ_exp.hg38.rda")


query.COAD_exp.hg38 <- GDCquery(project = "TCGA-COAD",
data.category = "Transcriptome Profiling",
data.type = "Gene Expression Quantification",
workflow.type = "HTSeq - FPKM-UQ",
barcode =  c("TCGA-14-0736-02A-01R-2005-01", "TCGA-06-0211-02A-02R-2005-01"))
GDCdownload(query.COAD_exp.hg38)
expdat <- GDCprepare(query = query.COAD_exp.hg38,
save = TRUE,
save.filename = "COAD_exp.hg38.rda")


#----------------------------------------#
# DNA methylation: aligned against hg38  #
#----------------------------------------#
query_met.hg38 <- GDCquery(project= "TCGA-COAD",
data.category = "DNA Methylation",
platform = "Illumina Human Methylation 450",
barcode = c("TCGA-AA-3712-01A-21D-1721-05","TCGA-CK-6747-01A-11D-1837-05","TCGA-AA-3502-11A-01D-1407-05","TCGA-D5-6536-01A-11D-1721-05","TCGA-CM-6676-01A-11D-1837-05",
"TCGA-G4-6306-01A-11D-1772-05","TCGA-AZ-6605-01A-11D-1837-05","TCGA-G4-6588-01A-11D-1772-05","TCGA-G4-6297-01A-11D-1721-05","TCGA-A6-2675-01A-02D-1721-05",
"TCGA-DM-A0X9-01A-11D-A153-05","TCGA-AA-3697-01A-01D-1721-05","TCGA-G4-6315-01A-11D-1721-05","TCGA-F4-6808-01A-11D-1837-05","TCGA-CM-5864-01A-01D-1651-05",
"TCGA-D5-6531-01A-11D-1721-05","TCGA-AD-6888-01A-11D-1926-05","TCGA-D5-6920-01A-11D-1926-05","TCGA-CA-6716-01A-11D-1837-05","TCGA-AZ-4614-01A-01D-1407-05",
"TCGA-5M-AATE-01A-11D-A40X-05","TCGA-CK-6748-01A-11D-1837-05","TCGA-D5-6537-01A-11D-1721-05","TCGA-AA-3506-11A-01D-1407-05","TCGA-3L-AA1B-01A-11D-A36Y-05",
"TCGA-AA-3495-11A-01D-1407-05","TCGA-G4-6320-01A-11D-1721-05","TCGA-AZ-6603-01A-11D-1837-05","TCGA-AZ-6601-11A-01D-1772-05","TCGA-D5-6538-01A-11D-1721-05",
"TCGA-DM-A282-01A-12D-A16X-05","TCGA-A6-5659-01A-01D-1651-05","TCGA-CK-6746-01A-11D-1837-05","TCGA-AA-3509-11A-01D-1407-05","TCGA-A6-2677-01A-01D-A27A-05"))
GDCdownload(query_met.hg38)
data.hg38 <- GDCprepare(query_met.hg38)


#----------------------------------------#
# DNA methylation: aligned against hg38  #
#----------------------------------------#
query_READ.met.hg38 <- GDCquery(project= "TCGA-READ",
data.category = "DNA Methylation",
platform = "Illumina Human Methylation 450",
barcode = c("TCGA-F5-6702-01A-11D-1828-05","TCGA-G5-6235-01A-11D-1734-05","TCGA-AG-A01Y-01A-41D-A081-05","TCGA-EI-7004-01A-11D-1926-05","TCGA-EI-6514-01A-11D-1734-05",
"TCGA-F5-6864-01A-11D-1926-05","TCGA-AG-A036-11A-11D-A081-05","TCGA-DC-6681-01A-11D-1828-05","TCGA-G5-6233-01A-11D-1734-05","TCGA-CI-6624-01C-11D-1828-05",
"TCGA-AG-3592-01A-02D-1734-05","TCGA-AG-A01W-01A-21D-A081-05","TCGA-CL-4957-01A-01D-1734-05","TCGA-AF-6672-01A-11D-1828-05","TCGA-CI-6619-01B-11D-1828-05",
"TCGA-AF-2690-01A-02D-1734-05","TCGA-F5-6465-01A-11D-1734-05","TCGA-AF-4110-01A-02D-1734-05","TCGA-EI-6507-01A-11D-1734-05","TCGA-DC-4749-01A-01D-1734-05",
"TCGA-AG-4021-01A-01D-1734-05","TCGA-EI-6884-01A-11D-1926-05","TCGA-EF-5831-01A-01D-1658-05","TCGA-DY-A0XA-01A-11D-A153-05","TCGA-G5-6572-02A-12D-1828-05",
"TCGA-F5-6863-01A-11D-1926-05","TCGA-BM-6198-01A-11D-1734-05","TCGA-F5-6814-01A-31D-1926-05","TCGA-DC-6683-01A-11D-1828-05","TCGA-AG-A01W-11A-11D-A081-05",
"TCGA-AF-A56N-01A-12D-A39G-05","TCGA-CI-6623-01B-11D-1828-05","TCGA-DC-5869-01A-01D-1658-05","TCGA-CI-6622-01A-11D-1828-05","TCGA-AG-3591-01A-01D-1734-05"))
GDCdownload(query_READ.met.hg38)
data.hg38 <- GDCprepare(query_READ.met.hg38)


#############################
## CNV data pre-processing ##
#############################
query.read_nocnv <- GDCquery(project = "TCGA-READ",
data.category = "Copy number variation",
legacy = TRUE,
file.type = "nocnv_hg19.seg",
sample.type = c("Primary solid Tumor"))
query.read_nocnv$results[[1]] <- query.read_nocnv$results[[1]][1:40,]

GDCdownload(query.read_nocnv, files.per.chunk = 100)

read_nocnv <- GDCprepare(query.read_nocnv, save = TRUE, save.filename = "READnocnvhg19.rda")


library(TCGAbiolinks)
# Obs: The data in the legacy database has been aligned to hg19
query.read_met <- GDCquery(project = "TCGA-READ",
legacy = TRUE,
data.category = "DNA methylation",
platform = "Illumina Human Methylation 450",
barcode = c("TCGA-F5-6702-01A-11D-1828-05","TCGA-G5-6235-01A-11D-1734-05","TCGA-AG-A01Y-01A-41D-A081-05","TCGA-EI-7004-01A-11D-1926-05","TCGA-EI-6514-01A-11D-1734-05",
"TCGA-F5-6864-01A-11D-1926-05","TCGA-AG-A036-11A-11D-A081-05","TCGA-DC-6681-01A-11D-1828-05","TCGA-G5-6233-01A-11D-1734-05","TCGA-CI-6624-01C-11D-1828-05",
"TCGA-AG-3592-01A-02D-1734-05","TCGA-AG-A01W-01A-21D-A081-05","TCGA-CL-4957-01A-01D-1734-05","TCGA-AF-6672-01A-11D-1828-05","TCGA-CI-6619-01B-11D-1828-05",
"TCGA-AF-2690-01A-02D-1734-05","TCGA-F5-6465-01A-11D-1734-05","TCGA-AF-4110-01A-02D-1734-05","TCGA-EI-6507-01A-11D-1734-05","TCGA-DC-4749-01A-01D-1734-05",
"TCGA-AG-4021-01A-01D-1734-05","TCGA-EI-6884-01A-11D-1926-05","TCGA-EF-5831-01A-01D-1658-05","TCGA-DY-A0XA-01A-11D-A153-05","TCGA-G5-6572-02A-12D-1828-05",
"TCGA-F5-6863-01A-11D-1926-05","TCGA-BM-6198-01A-11D-1734-05","TCGA-F5-6814-01A-31D-1926-05","TCGA-DC-6683-01A-11D-1828-05","TCGA-AG-A01W-11A-11D-A081-05",
"TCGA-AF-A56N-01A-12D-A39G-05","TCGA-CI-6623-01B-11D-1828-05","TCGA-DC-5869-01A-01D-1658-05","TCGA-CI-6622-01A-11D-1828-05","TCGA-AG-3591-01A-01D-1734-05"))
GDCdownload(query.read_met)

met.read_450 <- GDCprepare(query = query.read_met,
save = TRUE,
save.filename = "read_DNAmet450k.rda",
summarizedExperiment = TRUE)

query.met.coad <- GDCquery(project = "TCGA-COAD",
legacy = TRUE,
data.category = "DNA methylation",
platform = "Illumina Human Methylation 450",
barcode = c("TCGA-AA-3712-01A-21D-1721-05","TCGA-CK-6747-01A-11D-1837-05","TCGA-AA-3502-11A-01D-1407-05","TCGA-D5-6536-01A-11D-1721-05","TCGA-CM-6676-01A-11D-1837-05",
"TCGA-G4-6306-01A-11D-1772-05","TCGA-AZ-6605-01A-11D-1837-05","TCGA-G4-6588-01A-11D-1772-05","TCGA-G4-6297-01A-11D-1721-05",
"TCGA-DM-A0X9-01A-11D-A153-05","TCGA-AA-3697-01A-01D-1721-05","TCGA-G4-6315-01A-11D-1721-05","TCGA-F4-6808-01A-11D-1837-05","TCGA-CM-5864-01A-01D-1651-05",
"TCGA-D5-6531-01A-11D-1721-05","TCGA-AD-6888-01A-11D-1926-05","TCGA-D5-6920-01A-11D-1926-05","TCGA-CA-6716-01A-11D-1837-05","TCGA-AZ-4614-01A-01D-1407-05",
"TCGA-5M-AATE-01A-11D-A40X-05","TCGA-CK-6748-01A-11D-1837-05","TCGA-D5-6537-01A-11D-1721-05","TCGA-AA-3506-11A-01D-1407-05","TCGA-3L-AA1B-01A-11D-A36Y-05",
"TCGA-AA-3495-11A-01D-1407-05","TCGA-G4-6320-01A-11D-1721-05","TCGA-AZ-6603-01A-11D-1837-05","TCGA-AZ-6601-11A-01D-1772-05","TCGA-D5-6538-01A-11D-1721-05",
"TCGA-DM-A282-01A-12D-A16X-05","TCGA-CK-6746-01A-11D-1837-05","TCGA-AA-3509-11A-01D-1407-05"))
GDCdownload(query.met.coad)
met.coad.450 <- GDCprepare(query = query.met.coad,
save = TRUE,
save.filename = "coadDNAmet450k.rda",
summarizedExperiment = TRUE)
met.read_coad <- SummarizedExperiment::cbind(met.read_450, met.coad.450)


query.read_exp <- GDCquery(project = "TCGA-READ",
legacy = TRUE,
data.category = "Gene expression",
data.type = "Gene expression quantification",
platform = "Illumina HiSeq",
file.type = "results",
sample.type = "Primary solid Tumor")
GDCdownload(query.read_exp)
read_exp <- GDCprepare(query = query.read_exp, save = TRUE, save.filename = "read_Exp.rda")

query.exp.coad <- GDCquery(project = "TCGA-COAD",
legacy = TRUE,
data.category = "Gene expression",
data.type = "Gene expression quantification",
platform = "Illumina HiSeq",
file.type = "results",
sample.type = "Primary solid Tumor")
GDCdownload(query.exp.coad)
exp.coad <- GDCprepare(query = query.exp.coad, save = TRUE, save.filename = "exp.coad.rda")
exp.read_coad <- SummarizedExperiment::cbind(read_exp, exp.coad)

#-----------------------------------------------------------------------------
#                   Data.category: Copy number variation aligned to hg38
#-----------------------------------------------------------------------------
query.cnv <- GDCquery(project = "TCGA-READ",
data.category = "Copy Number Variation",
data.type = "Copy Number Segment",
barcode = c( "TCGA-OR-A5KU-01A-11D-A29H-01", "TCGA-OR-A5JK-01A-11D-A29H-01"))
GDCdownload(query)
data <- GDCprepare(query)

query <- GDCquery("TCGA-COAD",
"Copy Number Variation",
data.type = "Masked Copy Number Segment",
sample.type = c("Primary solid Tumor")) # see the barcodes with getResults(query)$cases
GDCdownload(query)
data <- GDCprepare(query)



#-----------------------------------------------------------------------------
#                    SummarizedExperiment object
#-----------------------------------------------------------------------------
library(SummarizedExperiment)

# Load object from TCGAWorkflowData package
# THis object will be created in the further sections,
data(READIllumina_HiSeq)

# get expression matrix
data <- assay(read_exp)
datatable(data[1:10,],
options = list(scrollX = TRUE, keys = TRUE, pageLength = 5),
rownames = TRUE)


#----------------------------
# Obtaining DNA methylation
#----------------------------
# Samples
read_samples <- matchedMetExp("TCGA-READ", n = 10)
coad.samples <- matchedMetExp("TCGA-COAD", n = 10)
samples <- c(read_samples,coad.samples)

#-----------------------------------
# 1 - Methylation
# ----------------------------------
# For methylation it is quicker in this case to download the tar.gz file
# and get the samples we want instead of downloading files by files
query <- GDCquery(project = c("TCGA-READ","TCGA-COAD"),
data.category = "DNA methylation",
platform = "Illumina Human Methylation 450",
legacy = TRUE,
barcode = samples)
GDCdownload(query)
met <- GDCprepare(query, save = FALSE)

# We will use only chr9 to make the example faster
met <- subset(met,subset = as.character(seqnames(met)) %in% c("chr9"))
# This data is avaliable in the package (object elmerExample)

#--------------------------
# DNA Methylation heatmap
#-------------------------
library(ComplexHeatmap)
clinical <- plyr::rbind.fill(read_clin,ccoad_clin)

# get the probes that are Hypermethylated or Hypomethylated
# met is the same object of the section 'DNA methylation analysis'
status.col <- "status.TCGA,READ.TCGA.COAD"
sig.met <- met[values(met)[,status.col], c("Hypermethylated","Hypomethylated"),]


# top annotation, which sampples are LGG and GBM
# We will add clinical data as annotation of the samples
# we will sort the clinical data to have the same order of the DNA methylation matrix
clinical.order <- clinical[match(substr(colnames(sig.met),1,12),clinical$bcr_patient_barcode),]
ta = HeatmapAnnotation(df = clinical.order[,c("disease","gender","vital_status","race")],
col = list(disease = c("READ" = "grey", "COAD" = "black"),
gender = c("male" = "blue","female" = "pink")))

# row annotation: add the status for LGG in relation to GBM
# For exmaple: status.gbm.lgg Hypomethyated means that the
# mean DNA methylation of probes for lgg are hypomethylated
# compared to GBM ones.
ra = rowAnnotation(df = values(sig.met)[status.col],
col = list("status.TCGA.READ.TCGA.COAD" =
c("Hypomethylated" = "orange",
"Hypermethylated" = "darkgreen")),
width = unit(1, "cm"))

heatmap  <- Heatmap(assay(sig.met),
name = "DNA methylation",
col = matlab::jet.colors(200),
show_row_names = FALSE,
cluster_rows = TRUE,
cluster_columns = FALSE,
show_column_names = FALSE,
bottom_annotation = ta,
column_title = "DNA Methylation")
# Save to pdf
png("heatmap.png",width = 600, height = 400)
draw(heatmap, annotation_legend_side =  "bottom")
dev.off()

#----------- 8.3 Identification of Regulatory Enhancers   -------
library(TCGAbiolinks)
# Samples: primary solid tumor w/ DNA methylation and gene expression
read_samples <- matchedMetExp("TCGA-READ", n = 10)
coad.samples <- matchedMetExp("TCGA-COAD", n = 10)
samples <- c(read_samples,coad.samples)

#-----------------------------------
# 1 - Methylation
# ----------------------------------
query.met <- GDCquery(project = c("TCGA-READ","TCGA-COAD"),
data.category = "DNA methylation",
platform = "Illumina Human Methylation 450",
legacy = TRUE,
barcode = samples)
GDCdownload(query.met)
met <- GDCprepare(query, save = FALSE)
met <- subset(met,subset = as.character(GenomicRanges::seqnames(met)) %in% c("chr9"))

#-----------------------------------
# 2 - Expression
# ----------------------------------
query.exp <- GDCquery(project = c("TCGA-LGG","TCGA-GBM"),
data.category = "Gene expression",
data.type = "Gene expression quantification",
platform = "Illumina HiSeq",
file.type  = "results",
legacy = TRUE,
barcode =  samples)
GDCdownload(query.exp)
exp <- GDCprepare(query.exp, save = FALSE)
save(exp, met, gbm.samples, lgg.samples, file = "elmer.example.rda", compress = "xz")

#----------------------------
# Obtaining DNA methylation
#----------------------------
# Samples
read_samples <- matchedMetExp("TCGA-READ", n = 10)
coad.samples <- matchedMetExp("TCGA-COAD", n = 10)
samples <- c(read_samples,coad.samples)


# We will use only chr9 to make the example faster
met <- subset(met,subset = as.character(seqnames(met)) %in% c("chr9"))
# This data is avaliable in the package (object elmerExample)
data(elmerExample)
#----------------------------
# Mean methylation
#----------------------------
# Plot a barplot for the groups in the disease column in the
# summarizedExperiment object

# remove probes with NA (similar to na.omit)
met <- subset(met,subset = (rowSums(is.na(assay(met))) == 0))

TCGAvisualize_meanMethylation(met,
groupCol = "project_id",
group.legend  = "Groups",
filename = NULL,
print.pvalue = TRUE)

met_DMR <- TCGAanalyze_DMR(met,
groupCol = "project_id", # a column in the colData matrix
group1 = "TCGA-READ", # a type of the disease type column
group2 = "TCGA-COAD", # a type of the disease column
p.cut = 0.05,
diffmean.cut = 0.15,
save = FALSE,
legend = "State",
plot.filename = "READ_COAD_metvolcano.png",
cores = 1)
