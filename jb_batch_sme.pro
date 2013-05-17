pro jb_batch_sme, catName, runNotes=runNotes
	; line-feed for pretty printing of multiple notes
	lf = string(10B)
	
	; the lockfile for the catalog is just the name of the catalog
	; with any extension removed and '.lock' appended
	dotPos = strpos(catName,'.',/REVERSE_SEARCH)
	if ( dotPos gt 0 ) then begin
		lockName = strmid(catName,0,dotPos) + ".lock"
	endif else lockName = catName + ".lock"
	
	; Get the process ID of this IDL session
	;	we will use this later to look for stop commands
	spawn,'/home/jbrewer/getppid.sh',pid
	stopFile = '/home/jbrewer/mdwarfs/iso_bat/'+pid+'.stop'
	
	if ( strpos(catName,'lick') ge 0 ) then begin
		savDir = '/home/jbrewer/mdwarfs/iso_lick/'
	endif else begin
		savDir = '/home/jbrewer/mdwarfs/iso_results/'
	endelse
	
	failedRun = 0
	
	repeat begin
		; Check to see if a note has been left for us to stop
		stopInfo = File_Info(stopFile)
		if ( stopInfo.Exists ) then exit
		
		; Check for existing catalog file and catalog lock
		;	wait for the lock to clear before opening catalog, take out a lock before editing.
		;
		ok = batchCatLock(lockName,runNotes)
		if ( ~ ok ) then stop, "Run Catalog Locked (" + catName + "), unable to get new runs."
		
		catInfo = File_Info(catName)
		if ( ~ catInfo.Exists ) then begin
			print, "Unable to locate run catalog: " + catName
			batchReleaseLock,lockName
			stop
		endif
		
		restore, catName ; jbIterCat
		unRunStars = where(jbIterCat.status eq 'NOT RUN',numAvailable)
		
		if (failedRun) then begin
			; The previous run by this processor crashed, so mark
			; that star as crashed before moving on.
			if ( n_elements(star) gt 0 ) then begin
				crashed = where(jbIterCat.obsnm eq star.obsnm AND jbIterCat.status eq 'RUNNING', nCrashed)
				if ( nCrashed eq 1 ) then begin
					jbitercat[crashed].status = 'CRASHED'
					jbitercat[crashed].runnotes = systime() + ": " + !ERROR_STATE.MSG + string(10B) + jbitercat[crashed].runnotes
					
					openw,ERRFILE,savDir+star.name+'_'+star.obsnm+'/error.txt',/GET_LUN
					printf,ERRFILE,systime() + ": " + !ERROR_STATE.MSG + string(10B) + jbitercat[crashed].runnotes
					free_lun,ERRFILE
				endif
			endif
			failedRun = 0
		endif
		
		;
		; If there are still stars to run, run the first one found.
		;
		if ( numAvailable gt 0 ) then begin
			print, "Found " + strtrim(string(numAvailable),2) + " stars to run."
			
			nextStar = unRunStars[0]
			star = jbIterCat[nextStar] ; run the next available star
			jbIterCat[nextStar].status = 'RUNNING'
			jbIterCat[nextStar].runStartDate = systime()
			if ( keyword_set(runNotes) ) then begin
				if ( jbIterCat[nextStar].runNotes ne '' ) then begin
					; add notes in reverse chronological order
					jbIterCat[nextStar].runNotes = runNotes + lf + jbIterCat[nextStar].runNotes
				endif else begin
					jbIterCat[nextStar].runNotes = runNotes
				endelse
			endif
			
			; Save the updated structure and release the lock
			save, jbIterCat, file=catName
			batchReleaseLock,lockName
			
			; Set an error handler so that if SME crashes, it doesn't take
			; the whole batch job down.
			;
			CATCH, Error_Status
			IF Error_status ne 0 then begin
				print, "ERROR running jb_sme_iter on " + star.name + " (" + star.obsnm + ") at " + systime()
				print, "ERRRO MSG: " + !ERROR_STATE.MSG
				failedRun = 1
				CATCH, /CANCEL
				continue ; skip to the next star, leaving this one at 'RUNNING'
				; The beginning of the loop will notice the error and set status to 'CRASHED'
			ENDIF
		
			; Run the star
			;
			print, "Starting SME Run for " + star.name + " (" + star.obsnm + ") at " + systime()
			if ( keyword_set(runNotes) ) then print, runNotes
			
			obNum = stregex(star.obsnm,'[^\.]+$',/EXTRACT)
			if ( strpos(catname,'lick') ge 0 ) then begin
				jb_sme_iter, star.runNum, star.name, obs=obNum, /lick
			endif else begin
				jb_sme_iter, star.runNum, star.name, obs=obNum
			endelse
			
			; Clear the error handler set above
			;
			CATCH, /CANCEL
			
			; Completed...we should check to see if this completed nicely before
			; just writing that it is done, but for now this will do.
			;
			ok = batchCatLock(lockName,"Finished Running " + star.name + " (" + star.obsnm + ")")
			if ( ~ ok ) then stop, "Unable to update catalog with completed star."
			
			; restore the catalog, find the star, update the status and run-end-time
			; save the catalog, clear the lock.
			;
			restore, catName
			thisStar = where(jbIterCat.obsnm eq star.obsnm and jbIterCat.status eq 'RUNNING',nFound)
			if ( nFound eq 0 ) then begin
				batchReleaseLock
				stop, "Running star not found, unable to update catalog."
			endif
			if ( nFound gt 1 ) then begin
				batchReleaseLock
				stop, "More than one running star found, unable to update catalog."
			endif
			jbIterCat[thisStar].status = 'COMPLETED'
			jbIterCat[thisStar].runEndDate = systime()
			save, jbIterCat, file=catName
			batchReleaseLock,lockName
		endif else begin
			; no more stars to run
			save, jbIterCat, file=catName ; in case we updated a crash status
			batchReleaseLock, lockName
		endelse
	endrep until ( numAvailable le 0 )
end