function parsed = irfImport(filename)
   fID=fopen(filename,'r','l');
   parsed=fread(fID,Inf,'double');
   fclose(fID);
end