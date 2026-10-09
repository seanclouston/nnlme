cap program drop nnlme
program define nnlme, eclass byable(onecall) prop(me xt mi bayes xtbs)
    version 17.0

    // Step 1: Parse syntax and options
    syntax varlist(min=1 max=1) [if] [in], ///  
        SURV(varlist fv) Time(varname)  ///
        [ID(varname)] ///
        [SHAPE(varlist fv)] [ITer(integer 25)] ///
        [   METHod(string)      ///
            RESCOVariance(string)   ///
            RESCORRelation(string)  ///
            RESVARiance(string) ///
            INITial(string)     ///
            lmeopts(string)     ///
            pnlsopts(string)    ///
            NOCONStant      ///
            xb          ///
            EMITERate(passthru)  ///
            EMTOLerance(passthru)   ///
            emlog           ///
        ]

    // Extract dependent variable
    gettoken depvar rest : varlist

    // Sample markout
    marksample touse

    // Validate required options
    if "`id'" == "" {
        di as error "option id() is required"
        exit 198
    }
    if "`surv'" == "" {
        di as error "option surv() is required"
        exit 198
    }
    if "`time'" == "" {
        di as error "option time() is required"
        exit 198
    }

    confirm numeric variable `time'

    // Display status
    di as txt " "
    di as txt ">>> Starting up NNLME <<<"
    di as txt " "
    
    if `iter' != 25 {
        di as txt ">>> Max Iteration Set: `iter'"
    }

    // Run MENL model wrapped in capture
    if "`shape'" != "" {
        capture noisily menl `depvar' = {ri1:} + {shp:} + {sl1:}*{surv:} if `touse', ///
            define(ri1: U0[`id'] , xb) ///
            define(shp: `shape' , xb nocons) ///
            define(sl1: -exp({a00})) ///
            define(surv: ln(1 + exp(min(700, `time' - exp(-{curv_c1: `surv', xb}))))) ///
            iterate(`iter') ///
            stddev noretable nofetable noheader pnlsopts(iterate(50) tolerance(1e-8)) emiterate(50)
    }
    else {
        capture noisily menl `depvar' = {ri1:} + {sl1:}*{surv:} if `touse', ///
            define(ri1: U0[`id'], xb) ///
            define(sl1: -exp({a00})) ///
            define(surv: ln(1 + exp(min(700, `time' - exp(-{curv_c1: `surv', xb}))))) ///
            iterate(`iter') ///
            stddev noretable nofetable noheader pnlsopts(iterate(50) tolerance(1e-8)) emiterate(50)
    }
        
    // Check execution
    local menl_rc = _rc
    if `menl_rc' {
        di as error ">>> Model estimation failed, code: `menl_rc'."
        exit `menl_rc'
    }

    // Extract e-returns safely using temporary matrix names
	matrix def BT = e(b)
    matrix def VT = e(V)
    
    local N = e(N)
    local L = e(ll)
    local fit = round(e(ll), .01)
    local p = round(e(p), .0001)
    local ch = round(e(chi2), .01)

    local group = .
    local hier = .
    capture mat H = e(hierstats)
    if !_rc {
        local group = H[1,1]
        local hier = round(H[1,2], .001)
    }

    local n_coef = colsof(BT)
    local eqs : coleq BT
    local names : colnames BT
    
    // Post results
    ereturn post BT VT, obs(`N') esample(`touse')
    ereturn scalar ll = `L'
    ereturn local cmd "nnlme"
   
    // Header Display
    di as txt _n "Nested Non-Linear Longitudinal Mixed-Effects Model" _col(58) "Number of obs" _col(75) "= " as res `N'
    di as txt "Group Variable: `id'" _col(58) "N. Groups" _col(75) "= " as res `group'
    di as txt " " _col(58) "Obs in groups" _col(75) "= " as res `hier'
    di as txt _col(58) "Wald chi2" _col(75) "= " as res "`ch'"
    di as txt "Linearization log likelihood = `fit'" _col(58) "Prob. > Chi2" _col(75) "= " as res "`p'"
    di as txt " " 
    
    // Tracking flags
    local surv_printed = 0
    local shape_printed = 0
    local a00_printed = 0
    local re_printed = 0
    local ri_printed = 0 
    local resid_printed = 0

    forvalues i = 1/`n_coef' {
        local eq : word `i' of `eqs'
        local vname : word `i' of `names'
		if strpos("`vname'", "b.") > 0 continue
		if strpos("`vname'", "o.") > 0 continue
        local display_vname "`vname'"

        // Identify parameter types
        local is_shape = ("`eq'" == "shp")
        local is_surv = ("`eq'" == "curv_c1")
        local is_a00 = ("`eq'" == "a00" | "`vname'" == "a00")
        local is_ri = ("`eq'" == "ri1")
        local is_resid = strpos("`vname'", "lnsigma")
        local is_re = (strpos("`eq'", "lnsd") | strpos("`vname'", "lnsd"))

        if !(`is_surv' | `is_shape' | `is_a00' | `is_re' | `is_resid' | `is_ri') continue

        // Print section headers
        if `is_surv' & !`surv_printed' {
            di as txt "{hline 80}"
            di as txt _col(14) "{c |}" _col(20) "NIR" _col(30) "Std. Err." _col(45) "Z" _col(52) "P>|Z|" _col(60) "[95% Conf. Interval]"
            di as txt "{hline 13}{c +}{hline 67}"    
            di as txt "/survival" _col(14) "{c |}"
            local surv_printed = 1
        }
        if `is_shape' & !`shape_printed' {
            di as txt "{hline 13}{c +}{hline 67}"
            di as txt _col(14) "{c |}" _col(20) "Coef." _col(30) "Std. Err." _col(45) "Z" _col(52) "P>|z|" _col(60) "[95% Conf. Interval]"
            di as txt "{hline 13}{c +}{hline 67}"
            di as txt "/shape" _col(14) "{c |}"
            local shape_printed = 1
        }
        if `is_a00' & !`a00_printed' {
            di as txt "{hline 13}{c +}{hline 67}"
            di as txt "Acc. param." _col(14) "{c |}"
            local a00_printed = 1
        }
        if `is_ri' & !`ri_printed' {
            di as txt "{hline 13}{c +}{hline 67}"
            di as txt "Intercept" _col(14) "{c |}"
            local ri_printed = 1
        }
        if `is_re' & !`re_printed' {
            di as txt "{hline 13}{c +}{hline 67}"
            di as txt "Random Eff." _col(14) "{c |}" _col(20) "Est." _col(30) "Std. Err." _col(45) "Z" _col(52) "P>|z|" _col(60) "[95% Conf. Interval]"
            di as txt "{hline 13}{c +}{hline 67}"
            di as txt "/`id'" _col(14) "{c |}"
            local re_printed = 1
        }
        if `is_resid' & !`resid_printed' {
            di as txt "/residual" _col(14) "{c |}"
            local resid_printed = 1
        }

        // Stats calculation using matrix element function
        tempname b_i v_i
        scalar `b_i' = el(e(b), 1, `i')
        scalar `v_i' = el(e(V), `i', `i')

        local coef   = `b_i'
        local se_i   = sqrt(`v_i')
        local z_i    = `coef' / `se_i'
        local p_i    = 2 * (1 - normal(abs(`z_i')))
        local ci_lb  = `coef' - 1.96 * `se_i'
        local ci_ub  = `coef' + 1.96 * `se_i'

        if `is_surv' {
            local display_coef  = exp(`coef')
            local display_se    = `display_coef' * `se_i'
            local display_ci_lb = exp(`ci_lb')
            local display_ci_ub = exp(`ci_ub')
        }   
        else {
            local display_coef  = `coef'
            local display_se    = `se_i'
            local display_ci_lb = `ci_lb'
            local display_ci_ub = `ci_ub'
        }

        // Format pre-calculated numbers directly into string macros
        local fmt_coef : display %9.0g `display_coef'
        local fmt_se   : display %9.0g `display_se'
        local fmt_z    : display %6.2f `z_i'
        local fmt_p    : display %6.3f `p_i'
        local fmt_lb   : display %9.0g `display_ci_lb'
        local fmt_ub   : display %9.0g `display_ci_ub'

        // Display pre-formatted text macros safely
        di as txt %12s "`vname'" _col(14) "{c |}" as res ///
            _col(16) "`fmt_coef'" ///
            _col(28) "`fmt_se'" ///
            _col(40) "`fmt_z'" ///
            _col(49) "`fmt_p'" ///
            _col(58) "`fmt_lb'" ///
            _col(70) "`fmt_ub'"
    }

    if e(converged) == 0 {
        di as txt "Warning: Convergence not achieved."
    }
    di as txt "{hline 80}"
end
