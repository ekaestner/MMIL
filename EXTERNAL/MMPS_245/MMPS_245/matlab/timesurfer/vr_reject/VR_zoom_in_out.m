function VR_zoom_in_out(varargin) %0 is zoom on main plot. 1 is zoom in raw. 2 is zoom on summary. 3 is zoom on raw via text
global CHANNELDATA
global TRIALDATA

src = varargin{1};
evnt= varargin{2};
code = varargin{3};
if code==0
    if evnt.VerticalScrollCount < 0
        child = get(CHANNELDATA.plotfig,'Children');
        panel = findobj(child,'flat','Type','uipanel','Visible','on');
        axes = get(panel,'Children');
        Ylim = get(axes(1),'Ylim');
        set(axes,'Ylim',Ylim/2);
        
    elseif evnt.VerticalScrollCount > 0
        child = get(CHANNELDATA.plotfig,'Children');
        panel = findobj(child,'flat','Type','uipanel','Visible','on');
        axes = get(panel,'Children');
        Ylim = get(axes(1),'Ylim');
        set(axes,'Ylim',Ylim*2);
        
    end;
elseif code==1 || code ==3;
    if code==1
        if evnt.VerticalScrollCount < 0
            child_raw = get(CHANNELDATA.raw_plot,'Children');
            panel_raw = findobj(child_raw,'flat','Type','uipanel','Visible','on');
            axes_raw = get(panel_raw,'Children');
            Ylim_raw = get(axes_raw(1),'Ylim');
            set(axes_raw,'Ylim',Ylim_raw/2);
        elseif evnt.VerticalScrollCount > 0
            child_raw = get(CHANNELDATA.raw_plot,'Children');
            panel_raw = findobj(child_raw,'flat','Type','uipanel','Visible','on');
            axes_raw = get(panel_raw,'Children');
            Ylim_raw = get(axes_raw(1),'Ylim');
            set(axes_raw,'Ylim',Ylim_raw*2);
        end;
    elseif code ==3
            Y_min = str2double(get(CHANNELDATA.raw_min,'String'));
            Y_max = str2double(get(CHANNELDATA.raw_max,'String'));
            child_raw = get(CHANNELDATA.raw_plot,'Children');
            panel_raw = findobj(child_raw,'flat','Type','uipanel','Visible','on');
            axes_raw = get(panel_raw,'Children');
            Ylim_raw = [Y_min Y_max];
            set(axes_raw,'Ylim',Ylim_raw); 
    end;
    Ylimit = get(axes_raw(1),'Ylim');
    set(CHANNELDATA.raw_min,'String',Ylimit(1));
    set(CHANNELDATA.raw_max,'String',Ylimit(2));

elseif code==2;
    panels=findobj(get(CHANNELDATA.sumplot,'Children'),'flat','Type','uipanel');
    if strcmp(get(panels(1),'Tag'),'trial_plots')
        t_axes=findobj(get(panels(1),'Children'),'Type','axes');
        c_axes=findobj(get(panels(2),'Children'),'Type','axes');
    else
        t_axes=findobj(get(panels(2),'Children'),'Type','axes');
        c_axes=findobj(get(panels(1),'Children'),'Type','axes');
    end;
    tlim = get(t_axes,'Ylim');
    clim = get(c_axes,'Xlim');
    
    if evnt.VerticalScrollCount < 0
        set(t_axes,'Ylim',tlim/2);
        set(c_axes,'Xlim',clim/2);
    elseif evnt.VerticalScrollCount > 0
        set(t_axes,'Ylim',tlim*2);
        set(c_axes,'Xlim',clim*2);
    end;
end;
end