dotprodcosine = function(df, query_mz_colname = 'mz', ref_mz_colname = 'theoretical.mz', 
                         query_int_colname = 'query_abund', ref_int_colname = 'ref_abund'){
  
  query_mz = df[,query_mz_colname]
  query_int = df[,query_int_colname]
  weighted_int_query = (query_mz^2) * sqrt(query_int)
  
  ref_mz = df[,ref_mz_colname]
  ref_int = df[,ref_int_colname]
  weighted_int_ref = (ref_mz^2) * sqrt(ref_int)
  
  numerator = sum(weighted_int_query * weighted_int_ref)
  denominator = sqrt( sum(weighted_int_query^2) ) * sqrt( sum(weighted_int_ref^2) )
  
  dotprodcosine = numerator / denominator
  
  return(dotprodcosine)
  
}