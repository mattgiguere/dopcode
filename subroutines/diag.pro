pro diag, units, text, prefix=prefix, wrap=wrap, stop=stop, retall=retall
;
;Print a multi-line diagnostic message to the specifed file units.
;
;Mandatory Inputs:
;
; units (scalar or vector integer) list of logical file units where text
;   will be written. unit -1 is standard output. unit -2 is standard error.
;   unit -3 suppresses output.
;
; text (scalar or vector string) lines of text that will be written to units.
;
;Optional Inputs:
;
; [prefix=] (string) prefix to write before each line of text. if prefix is
;   not specifed, then is set to '<routine_name>@<line_number>: '
;
; [wrap=] (integer) maximum number of characters to output per output line
;   before wrapping to a new line. line breaks will occur at word boudaries,
;   if possible. word boundaries are defined as spaces in the output string.
;
; [/stop] (switch) causes program control to halt in the calling routine
;   after witing the specified text.
;
; [/retall] (switch) causes program control to return to the top level after
;   the specified text is written.
;
;History:
;
; 2008-Feb-29 Valenti  Coded initial version.
; 2008-Dec-06 Valenti  Added line wrap logic.

;Check syntax.
  if n_params() lt 2 then begin
    print, 'syntax: diag, units, text [,prefix= ,wrap= ,/stop ,/retall]'
    print, "  e.g.: diag, [-1,unit], pre='', 'file not found', /ret"
    return
  endif

;Default values for optional parameters.
  if not keyword_set(wrap) then wrap = 999

;Verify that unit is specified.
  nunit = n_elements(units)
  if nunit eq 0 then message, 'output unit(s) not defined'

;Exit without printing output, if units has flag value of -3.
  if nunit eq 1 and units[0] eq -3 then return

;Get traceback information, if needed.
  if n_elements(prefix) eq 0 or keyword_set(return) then begin
;
;   callstack = scope_traceback(/struct)		;does not compile...
;   caller = callstack[n_elements(callstack)-2]		;...for version < 6.2
;
    help, /trace, output=help_trace
    words = strsplit(help_trace[1], /extract)
    if n_elements(words) eq 2 then begin
      caller = { routine: words[1], line: -1 }
    endif else begin
      caller = { routine: words[1], line: long(words[2]) }
    endelse
  endif

;If prefix is not specified, set prefix to the name of the calling routine.
  if n_elements(prefix) eq 0 then begin
    if caller.routine eq '$MAIN$' then begin
      prefix = 'main'
    endif else begin
      prefix = strlowcase(caller.routine) + '@' + strtrim(caller.line, 2)
    endelse
    prefix = prefix + ': '
  endif

;Loop through lines of output text, printing to every specified unit.
;Wrap long lines, preferably at word boundaries.
  maxlen = wrap - strlen(prefix)
  for itext=0, n_elements(text)-1 do begin		;loop thru output lines
    out = text[itext]					;extract current line
    while strlen(out) gt maxlen do begin		;true: need to wrap
      ipos = strpos(out, ' ', maxlen, /reverse_search)	; find word boundary
      if ipos eq -1 then ipos = maxlen			; no space, so cut word
      for iunit=0, nunit-1 do begin			; loop thru units
        printf, units[iunit], prefix + strmid(out, 0, ipos)	;print portion
      endfor
      out = strmid(out, ipos, 999999)			; trim printed portion
    endwhile
    for iunit=0, nunit-1 do printf, units[iunit], prefix + out	;print remains
  endfor

;Stop execution in calling routine, if requested.
  if keyword_set(stop) then begin
    on_error, 2					;return to caller on error
    message, /noname, /noprefix, ''
  endif

;Break to top level context, if requested.
  if keyword_set(retall) then retall

end
