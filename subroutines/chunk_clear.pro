pro chunk_clear, chunk_arr

nel = n_elements(chunk_arr)
tags = tag_names(chunk_arr)
ntags = n_elements(tags)

for i=0,nel - 1 do begin
   for j = 0, ntags-1 do begin 
      if size(chunk_arr[i].(j), /type) eq 10 then begin 
	ptr_free, chunk_arr[i].(j)
    endif
   endfor
endfor

end
