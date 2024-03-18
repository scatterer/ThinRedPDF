
function y_reduced = bkgr_reduction(x_scan, y_scan, x_bkgr, y_bkgr)

if 1==1
    %compensation_factor = y_scan(1)/y_bkgr(1);
    compensation_factor = 0.5; %fminsearch(@find_bkgr_compensation_factor, 0.5);
    y_reduced = y_scan - compensation_factor.*y_bkgr;
else
    y_reduced = 1;
end

    function I_red = find_bkgr_compensation_factor(bkgr_factor)
        index = 10;
        I_red =  std( y_scan(1:index) -  bkgr_factor.*y_bkgr(1:index) );
        
    end

end