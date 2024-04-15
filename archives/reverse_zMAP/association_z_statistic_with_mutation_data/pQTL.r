library(MatrixEQTL)
useModel = modelLINEAR;

output_file_name = tempfile();

pvOutputThreshold = 0.1;

errorCovariance = numeric();


snps = SlicedData$new();
snps$fileDelimiter = "\t";
snps$fileOmitCharacters = "NA";
snps$fileSkipRows = 1;
snps$fileSkipColumns = 1;
snps$fileSliceSize = 2000;
snps$LoadFile("mutation_df.txt");


gene = SlicedData$new();
     gene$fileDelimiter = "\t";      # the TAB character
     gene$fileOmitCharacters = "NA"; # denote missing values;
     gene$fileSkipRows = 1;          # one row of column labels
     gene$fileSkipColumns = 1;       # one column of row labels
     gene$fileSliceSize = 2000;      # read file in slices of 2,000 rows
gene$LoadFile("z_df.txt");


cvrt = SlicedData$new();
 cvrt$fileDelimiter = "\t";      # the TAB character
 cvrt$fileOmitCharacters = "NA"; # denote missing values;
 cvrt$fileSkipRows = 1;          # one row of column labels
 cvrt$fileSkipColumns = 1;       # one column of row labels
 if(length("covariates.txt")>0) {
 cvrt$LoadFile("covariates.txt");
 }


me = Matrix_eQTL_engine(
     snps = snps,
   gene = gene,
   cvrt = cvrt,
   output_file_name = output_file_name,
   pvOutputThreshold = pvOutputThreshold,
   useModel = useModel,
   errorCovariance = errorCovariance,
   verbose = TRUE,
   pvalue.hist = TRUE,
   min.pv.by.genesnp = FALSE,
   noFDRsaveMemory = FALSE);


unlink(output_file_name);


#cat('Analysis done in: ', me$time.in.sec, ' seconds', '\n');


#cat('Detected eQTLs:', '\n');


#show(me$all$eqtls)


plot(me)


results <- me$all$eqtls


write.table (results, file ="results.txt", sep ="\t",row.names =FALSE, col.names =TRUE,quote=FALSE)
