function parameters = initalise_parameters(ns_ps, Srate, method)

len_val = length(ns_ps)

switch lower(method)
  case 'mcra'
    parameters = struct('n',2,'len',len_val,'P',ns_ps,'Pmin',ns_ps,'Ptmp',ns_ps,'pk',zeros(len_val,1),'noise_ps',ns_ps,...
            'ad',0.95,'as',0.8,'L',round(1000*2/20),'delta',5,'ap',0.2);
end

end % function parameters = initalise_parameters(ns_ps, Srate, method)