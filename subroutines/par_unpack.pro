function par_unpack, chunk, names=names, ind=ind
struct = (*chunk.par)[0]
tags = tag_names(struct)
ntags = n_elements(tags)
ind = intarr(2, ntags)
eachname = ['WCOF', 'AMP', 'OFFSET', 'Z']
nen = n_elements(eachname)
if ntags gt nen then eachname = [eachname, strarr(ntags-nen)+'EXTRA']
for i = 0, ntags-1 do begin
    nel = n_elements(struct.(i))
    if i eq 0 then begin
        ind[*,i] = [0, n_elements(struct.(i))-1]
        par = struct.(i) 
        nums = strtrim(string(indgen(nel), format=format),2)
        names = eachname[i]+nums
    endif else begin
        ind[*,i] = ind[1,i-1] + [0, n_elements(struct.(i))-1]
        par = [par, struct.(i)]
        nums = strtrim(string(indgen(nel), format=format),2)
        names = [names, eachname[i]+nums]
    endelse
endfor
return, par
end
