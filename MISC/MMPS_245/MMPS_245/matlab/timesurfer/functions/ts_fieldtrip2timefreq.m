function data = ts_fieldtrip2timefreq(ft_data, template_data)
%function data = ts_fieldtrip2timefreq(ft_data, template_data)

  if (~exist('template_data','var'))
    data = ts_fieldtrip2data(ft_data, 'timefreq');
  else
    data = ts_fieldtrip2data(ft_data, 'timefreq', template_data);
  end;