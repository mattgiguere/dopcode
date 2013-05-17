function check_file,fname

ff=findfile(fname,count=exists)
if exists gt 0 then exists=1

return,exists
end
