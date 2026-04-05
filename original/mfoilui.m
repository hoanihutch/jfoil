function mfoilui
%-------------------------------------------------------------------------------
% mfoilui.m: user interface for mfoil class (v 2021-10-10)
%
% Copyright (C) 2021 Krzysztof J. Fidkowski
%
% Permission is hereby granted, free of charge, to any person
% obtaining a copy of this software and associated documentation files
% (the "Software"), to deal in the Software without restriction,
% including without limitation the rights to use, copy, modify, merge,
% publish, distribute, sublicense, and/or sell copies of the Software,
% and to permit persons to whom the Software is furnished to do so,
% subject to the following conditions:
%
% The above copyright notice and this permission notice shall be
% included in all copies or substantial portions of the Software.
%
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
% EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
% MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
% NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS
% BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN
% ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN
% CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%-------------------------------------------------------------------------------

% default mfoil object
snaca = '2412';
M = mfoil('naca',snaca); 
M.param.doplot = false;
  
% default row height
RH = 20;
  
% use a fixed-size figure
P = [100, 100, 900, 600];

% create UI figure
fig = uifigure('name',sprintf('mfoil (version %s)', M.version),'windowstyle','normal');
P = get(fig,'position'); P(3) = 900; P(4) = 600; set(fig, 'position',P);

% top-level layout
gl0 = uigridlayout(fig,[3 3]);
%gl0.RowHeight = {'2x','2x','1x'};
gl0.RowHeight = {'4x','3x','2x'};
gl0.ColumnWidth = {'2x','2x','5x'};

% gl0 components: panels and grid layouts withing panels
% airfoil
pan_foil  = uipanel(gl0);
gl_foil = uigridlayout(pan_foil, [6,4]);
gl_foil.RowHeight = num2cell(ones(1,6)*RH);
% oper
pan_oper  = uipanel(gl0);
gl_oper = uigridlayout(pan_oper, [7,4]);
gl_oper.RowHeight = num2cell(ones(1,8)*RH);
%gl_oper.ColumnWidth = {'1x','1x'};
% param
pan_param  = uipanel(gl0);
gl_param = uigridlayout(pan_param, [6,2]);
gl_param.RowHeight = num2cell(ones(1,6)*RH);
gl_param.ColumnWidth = {'2x','1x'};
% plot
pan_plot  = uipanel(gl0);
gl_plot = uigridlayout(pan_plot, [4, 2]);
gl_plot.RowHeight = {RH, RH, 10, 112};
gl_plot.RowSpacing = 5;
% axes
ax_plot = uiaxes(gl0); ax_foil = uiaxes(gl0); ax_pan  = uiaxes(gl0);
axs = {ax_plot, ax_foil, ax_pan};
for i = 1:length(axs), axs{i}.NextPlot = 'replace'; axis(axs{i},'off'); end;

% component layout
pan_foil.Layout.Row = 1;     pan_foil.Layout.Column = 1;
pan_oper.Layout.Row = 1;     pan_oper.Layout.Column = 2;
pan_param.Layout.Row = 2;    pan_param.Layout.Column = 1;
pan_plot.Layout.Row = 2;     pan_plot.Layout.Column = 2;
ax_plot.Layout.Row = [1,2];  ax_plot.Layout.Column = 3;
ax_foil.Layout.Row = 3;      ax_foil.Layout.Column = 3;
ax_pan.Layout.Row = 3;       ax_pan.Layout.Column = [1,2];

% all panels and axes
panels = {pan_foil, pan_oper, pan_param, pan_plot, ax_plot, ax_foil, ax_pan};


% -- airfoil panel --
% name
foil_name_lb = uilabel(gl_foil, 'Tooltip', 'Name of the current airfoil');
set(foil_name_lb, 'HorizontalAlignment','right', 'Text','Aifoil name');
foil_name_lb.Layout.Column = [1,2];
foil_name_ef = uieditfield(gl_foil,'Value', M.geom.name, ...
                           'ValueChangedFcn', @(txt,event) foil_name_set(M,txt));
foil_name_ef.Layout.Column = [3,4];
% naca
foil_naca_lb = uilabel(gl_foil, 'Tooltip', 'Load a 4 or 5-digit NACA airfoil');
set(foil_naca_lb, 'HorizontalAlignment', 'right', 'Text', 'NACA');
foil_naca_lb.Layout.Column = [1,2];
foil_naca_ef = uieditfield(gl_foil,'Value', snaca, ...
                           'ValueChangedFcn', @(txt,event) foil_naca_set(M,panels,txt));
foil_naca_ef.Layout.Column = [3,4];
% file load
foil_load_lb = uilabel(gl_foil, 'Tooltip', 'Load a plain-text coordinate file');
set(foil_load_lb, 'HorizontalAlignment', 'right', 'Text', 'Coordinate file');
foil_load_lb.Layout.Column = [1,2];
foil_load_bt = uibutton(gl_foil, 'Text', 'Load', 'Tooltip', 'Load a plain-text coordinate file', ...
                       'ButtonPushedFcn', @(btn,event) foil_load_press(M,panels,btn));
foil_load_bt.Layout.Column = [3,4];
% npanel
foil_npanel_lb = uilabel(gl_foil, 'Tooltip', 'Set number of panels on the airfoil');
set(foil_npanel_lb, 'HorizontalAlignment', 'right', 'Text', '# panels');
foil_npanel_lb.Layout.Column = [1,2];
foil_npanel_ef = uieditfield(gl_foil,'numeric', 'limits', [0,999], 'Value', M.foil.N-1, ...
                           'ValueChangedFcn', @(num,event) foil_npanel_set(M,panels,num));
foil_npanel_ef.Layout.Column = [3,4];
% flap
foil_flap_bt = uibutton(gl_foil, 'Text', 'Flap', 'Tooltip', 'Deploy a flap on the current airfoil', ...
                        'ButtonPushedFcn', @(btn,event) foil_flap_press(M,panels,btn));
flap_xhinge_ef = uieditfield(gl_foil,'numeric', 'Value', 0.8, 'Tooltip', 'flap x-hinge location', ...
                           'ValueChangedFcn', @(num,event) plot_panels(M,panels));
flap_zhinge_ef = uieditfield(gl_foil,'numeric', 'Value', 0.0, 'Tooltip', 'flap z-hinge location', ...
                           'ValueChangedFcn', @(num,event) plot_panels(M,panels));
% TODO: use inf for flap y value to mean middle
flap_angle_ef = uieditfield(gl_foil,'numeric', 'Value', 0.0, 'Tooltip', 'flap angle [deg]: positive is down');
% xref
foil_xref_lb = uilabel(gl_foil, 'Tooltip', 'Set moment reference coordinate');
set(foil_xref_lb, 'HorizontalAlignment', 'right', 'Text', 'Mom ref point');
foil_xref_lb.Layout.Column = [1,2];
xref_x_ef = uieditfield(gl_foil,'numeric', 'Value', M.geom.xref(1), 'Tooltip', 'moment center x location', ...
                        'ValueChangedFcn', @(num,event) set_xref_x(M,panels,num));
xref_z_ef = uieditfield(gl_foil,'numeric', 'Value', M.geom.xref(2), 'Tooltip', 'moment center z location', ...
                        'ValueChangedFcn', @(num,event) set_xref_z(M,panels,num));
% camber increment
foil_camb_lb = uilabel(gl_foil, 'Tooltip', 'Add a plain-text camber file');
set(foil_camb_lb, 'HorizontalAlignment', 'right', 'Text', 'Camber file');
foil_camb_lb.Layout.Column = [1,2];
foil_camb_bt = uibutton(gl_foil, 'Text', 'add camber', 'Tooltip', 'Load a plain-text file with camber addition', ...
                       'ButtonPushedFcn', @(btn,event) foil_camb_press(M,panels,btn));
foil_camb_bt.Layout.Column = [3,4];


%-- oper panel ---
% angle of attack
oper_alpha_lb = uilabel(gl_oper, 'Tooltip', 'Target angle of attack for alpha-constrained runs');
set(oper_alpha_lb, 'HorizontalAlignment', 'right', 'Text', '\alpha [deg]', 'Interpreter', 'tex');
oper_alpha_lb.Layout.Column = [1,2];
oper_alpha_ef = uieditfield(gl_oper,'numeric', 'limits', [-90 90], 'Value', M.oper.alpha, ...
                           'ValueChangedFcn', @(num,event) oper_alpha_set(M,num));
oper_alpha_ef.Layout.Column = [3,4];

% target lift coefficient
oper_cl_lb = uilabel(gl_oper, 'Tooltip', 'Target lift coefficient for lift-constrained runs');
set(oper_cl_lb, 'HorizontalAlignment', 'right', 'Text', 'target $c_{\ell}$', 'Interpreter', 'latex');
oper_cl_lb.Layout.Column = [1,2];
oper_cl_ef = uieditfield(gl_oper,'numeric', 'limits', [-10 10], 'Value', M.oper.cltgt, ...
                           'ValueChangedFcn', @(num,event) oper_cl_set(M,gl_oper,num));
oper_cl_ef.Layout.Column = [3,4];

% Reynolds number
oper_Re_lb = uilabel(gl_oper, 'Tooltip', 'Relevant for viscous flows only, based on the chord');
set(oper_Re_lb, 'HorizontalAlignment', 'right', 'Text', 'Reynolds #');
oper_Re_lb.Layout.Column = [1,2];
oper_Re_ef = uieditfield(gl_oper,'numeric', 'limits', [0 1e10], 'Value', M.oper.Re, ...
                           'ValueChangedFcn', @(num,event) oper_Re_set(M,gl_oper,num));
oper_Re_ef.Layout.Column = [3,4];

% Mach number
oper_Ma_lb = uilabel(gl_oper, 'Tooltip', 'Set > 0 to use the compressibility correction');
oper_Ma_lb.Layout.Column = [1,2];
set(oper_Ma_lb, 'HorizontalAlignment', 'right', 'Text', 'Mach #');
oper_Ma_ef = uieditfield(gl_oper,'numeric', 'limits', [0 1e10], 'Value', M.oper.Ma, ...
                           'ValueChangedFcn', @(num,event) oper_Ma_set(M,num));
oper_Ma_ef.Layout.Column = [3,4];


% Forced transition locations
oper_xftl_lb = uilabel(gl_oper, 'Tooltip', 'Forced transition x/c on lower surface');
oper_xftl_lb.Layout.Column = [1];
set(oper_xftl_lb, 'HorizontalAlignment', 'right', 'Text', 'xftl');
oper_xftl_ef = uieditfield(gl_oper,'numeric', 'limits', [0 1], 'Value', M.oper.xft(1), ...
                           'ValueChangedFcn', @(num,event) oper_xftl_set(M,num));
oper_xftl_ef.Layout.Column = [2];

oper_xftu_lb = uilabel(gl_oper, 'Tooltip', 'Forced transition x/c on upper surface');
oper_xftu_lb.Layout.Column = [3];
set(oper_xftu_lb, 'HorizontalAlignment', 'right', 'Text', 'xftu');
oper_xftu_ef = uieditfield(gl_oper,'numeric', 'limits', [0 1], 'Value', M.oper.xft(2), ...
                           'ValueChangedFcn', @(num,event) oper_xftu_set(M,num));
oper_xftu_ef.Layout.Column = [4];


% viscous state button
oper_viscous_bt = uibutton(gl_oper, 'state', 'Text', 'Inviscid', ...
                           'Tooltip', 'Toggle between inviscid/viscous modes', ...
                           'ValueChangedFcn', @(btn,event) oper_viscous_press(M,btn));
oper_viscous_bt.Layout.Column = [1,2];

% target cl button
oper_givencl_bt = uibutton(gl_oper, 'state', 'Text', 'Target alpha', ...
                           'Tooltip', 'Toggle between alpha/cl constrained modes', ...
                           'ValueChangedFcn', @(btn,event) oper_givencl_press(M,btn));
oper_givencl_bt.Layout.Column = [3,4];

% init BL button
oper_initbl_bt = uibutton(gl_oper, 'state', 'Text', 'Init BL', ...
                          'Tooltip', 'Toggle between initializing and reusing boundary-layer solution', ...
                          'ValueChangedFcn', @(btn,event) oper_initbl_press(M,btn));
oper_initbl_bt.Layout.Column = [1,2];

% run button
oper_run_bt = uibutton(gl_oper, 'Text', 'Run', 'FontWeight', 'bold', ...
                       'BackgroundColor', [0.2 0.9 0.2], 'Tooltip', 'Run the simulation', ...
                       'ButtonPushedFcn', @(btn,event) oper_run_press(M,panels,axs,btn));
oper_run_bt.Layout.Column = [3,4];

oper_conv_lb = uilabel(gl_oper, 'Tooltip', 'Indicates whether or not the coupled solver converged');
set(oper_conv_lb, 'HorizontalAlignment', 'center', 'Text', '');
oper_conv_lb.Layout.Column = [1,4];


%-- param panel --
% Newton iterations
param_niglob_lb = uilabel(gl_param, 'Tooltip', 'Number of coupled-solver Newton-Raphson iterations');
set(param_niglob_lb, 'HorizontalAlignment', 'right', 'Text', 'Newton iterations');
param_niglob_ef = uieditfield(gl_param,'numeric', 'limits', [1 500], 'Value', M.param.niglob, ...
                           'ValueChangedFcn', @(num,event) param_niglob_set(M,num));
% residual tolerance
param_rtol_lb = uilabel(gl_param, 'Tooltip', 'Residual termination criterion for the coupled solver');
set(param_rtol_lb, 'HorizontalAlignment', 'right', 'Text', 'Residual tolerance');
param_rtol_ef = uieditfield(gl_param,'numeric', 'limits', [0 2], 'Value', M.param.rtol, ...
                           'ValueChangedFcn', @(num,event) param_rtol_set(M,num));
% ncrit
param_ncrit_lb = uilabel(gl_param, 'Tooltip', 'This controls transition in the e^n model');
set(param_ncrit_lb, 'HorizontalAlignment', 'right', 'Text', 'Critical amplification');
param_ncrit_ef = uieditfield(gl_param,'numeric', 'limits', [0 20], 'Value', M.param.ncrit, ...
                           'ValueChangedFcn', @(num,event) param_ncrit_set(M,num));
% verbosity
param_verb_lb = uilabel(gl_param, 'Tooltip', 'Controls how much is printed to the terminal (4 = most)');
set(param_verb_lb, 'HorizontalAlignment', 'right', 'Text', 'Verbosity');
param_verb_ef = uieditfield(gl_param,'numeric', 'limits', [0 4], 'Value', M.param.verb, ...
                            'RoundFractionalValues', 'on', ...
                            'ValueChangedFcn', @(num,event) param_verb_set(M,num));

%-- plot panel --
% plotting drop-down
plot_dist_lb = uilabel(gl_plot, 'HorizontalAlignment', 'left', 'Text', 'Quantity to plot:');
plot_dist_lb.Layout.Column = [1,2];
dist_items = {'pressure coefficient', ...
              'skin friction coefficient', ...
              'displacement thickness', ...
              'momentum thickness', ...
              'kinematic shape parameter', ...
              'edge velocity', ...
              'amplification/shear-lag', ...
              'Reynolds_theta'};
plot_dist_dd = uidropdown(gl_plot, 'Items', dist_items, 'Value', 1, ...
                          'ItemsData', 1:length(dist_items), ...
                          'ValueChangedFcn', @(dd,event) plot_dist_set(M,axs,dd));
plot_dist_dd.Layout.Row = [2]; plot_dist_dd.Layout.Column = [1,2];
T = get_post_table(M);
plot_post_t = uitable(gl_plot, 'Data', T, 'ColumnName', [], 'ColumnWidth', {'2x','1x'}, ...
                      'Tooltip', 'Outputs calculated from the latest run');
plot_post_t.Layout.Row = [4]; plot_post_t.Layout.Column = [1,2];

% plot panels
plot_panels(M,panels);

end


%===============================================================================
% callbacks and supporting functions

%-------------------------------------------------------------------------------
function foil_name_set(M,txt)
  M.geom.name = txt.Value;
end

%-------------------------------------------------------------------------------
function foil_naca_set(M,panels,txt)
  M.geom_change('naca',txt.Value);
  pan_foil = panels{1};  gl_foil = pan_foil.Children(1); 
  foil_load_bt = gl_foil.Children(6); % load button
  foil_load_bt.Text = 'Load';
  clear_solution(M,panels);
end

%-------------------------------------------------------------------------------
function foil_npanel_set(M,panels,num)
  M.geom_panel(num.Value);
  clear_solution(M,panels);    
end

%-------------------------------------------------------------------------------
function foil_load_press(M,panels,btn)
  [file,path] = uigetfile('*.*', 'Choose a plain-text coordinate file');
  ffile = fullfile(path,file);
  X = load(ffile); % attempt to load text file
  M.geom_change('coords',X);
  M.geom.name = file; btn.Text = file;
  pan_foil = panels{1};  gl_foil = pan_foil.Children(1); 
  foil_naca_ef = gl_foil.Children(4); % naca field
  foil_naca_ef.Value = '';
  clear_solution(M,panels);
end

%-------------------------------------------------------------------------------
function foil_camb_press(M,panels,btn)
  [file,path] = uigetfile('*.*', 'Choose a plain-text camber file');
  ffile = fullfile(path,file);
  X = load(ffile); % attempt to load text file
  M.geom_addcamber(X);
  M.geom.name = [M.geom.name, ' camb'];
  clear_solution(M,panels);
end

%-------------------------------------------------------------------------------
function clear_solution(M,panels)
  pan_foil = panels{1};
  gl_foil = pan_foil.Children(1);
  foil_name_ef = gl_foil.Children(2); % name field
  foil_name_ef.Value = M.geom.name;
  pan_plot = panels{4};
  gl_plot = pan_plot.Children(1);
  gl_plot.Children(3).Data = get_post_table(M);
  cla(panels{end-2}); cla(panels{end-1}); cla(panels{end});
  plot_panels(M,panels);
end

%-------------------------------------------------------------------------------
function foil_flap_press(M,panels,btn)
  pan_foil = panels{1};  gl_foil = pan_foil.Children(1); 
  xzhinge = [gl_foil.Children(10).Value, gl_foil.Children(11).Value];
  eta = gl_foil.Children(12).Value;
  if (eta == 0), return; end; % nothing to do
  M.geom_flap(xzhinge, eta);
  M.geom.name = [M.geom.name, ' flap'];
  clear_solution(M,panels);
end

%-------------------------------------------------------------------------------
function set_xref_x(M,panels,num)
  M.geom.xref(1) = num.Value;
  clear_solution(M,panels);
end

%-------------------------------------------------------------------------------
function set_xref_z(M,panels,num)
  M.geom.xref(2) = num.Value;
  clear_solution(M,panels);
end

%-------------------------------------------------------------------------------
function plot_panels(M,panels)
  ax = panels{end};
  M.param.axplot = ax; cla(ax); ax.NextPlot='replace';
  chord = M.geom.chord;
  xz = [M.foil.x]; x = xz(1,:); N = M.foil.N;
  plot(ax, xz(1,:), xz(2,:), 'bo-', 'linewidth', 1, 'markersize', 3);
  xrange = [min(x)-0.01*chord, max(x)+0.01*chord]; 
  yrange = [min(M.foil.x(2,:)), max(M.foil.x(2,:))];
  V = [xrange, yrange];
  AR = 80/380; % panel aspect ratio
  dy = V(4)-V(3); dx = xrange(2)-xrange(1); V(3:4) = V(3:4)*(dx*AR)/dy;
  axis(ax,V); grid(ax, 'on'); %box(ax, 'off'); axis(ax,'off');
  pan_foil = panels{1};  gl_foil = pan_foil.Children(1); 
  hold(ax,'on'); 
  % flap hinge
  xzhinge = [gl_foil.Children(10).Value, gl_foil.Children(11).Value];
  hp(1) = plot(ax, xzhinge(1), xzhinge(2), 'rx', 'linewidth', 2, 'DisplayName', 'Flap hinge');
  % moment reference location
  xref = [gl_foil.Children(14).Value, gl_foil.Children(15).Value];
  hp(2) = plot(ax, xref(1), xref(2), 'm*', 'linewidth', 2, 'DisplayName', 'Moment point');
  % transition location
  if (M.oper.viscous)
    for is = 1:2
      xt = M.vsol.Xt(is,2);
      if (xt > 0)
        Is = M.vsol.Is{is};
        zt = interp1(xz(1,Is), xz(2,Is), xt, 'linear');
        hp(3) = plot(ax, xt, zt, 'g|', 'linewidth', 2, 'DisplayName', 'Transition');
      end
    end
  end
  h = legend(ax, hp, 'Location', 'North'); set(h, 'fontsize', 8, 'box', 'off', 'orientation', 'horizontal');
  %P = get(h,'Position'); P(1) = 0; P(2) = P(2)*1.1; P(3) = P(3)*0.8; set(h,'Position', P); get(h)
end

%-------------------------------------------------------------------------------
function oper_alpha_set(M,num)
  M.setoper('alpha',num.Value);
end

%-------------------------------------------------------------------------------
function oper_cl_set(M,gl_oper,num)
  M.setoper('cl',num.Value);
  h = gl_oper.Children(10);
  if (~h.Value), oper_givencl_set(M,h); end
end

%-------------------------------------------------------------------------------
function oper_Re_set(M,gl_oper,num)
  M.setoper('Re',num.Value);
  h = gl_oper.Children(9);
  if (~h.Value), oper_viscous_set(M,h); end
end

%-------------------------------------------------------------------------------
function oper_Ma_set(M,num)
  M.setoper('Ma',num.Value);
end

%-------------------------------------------------------------------------------
function oper_xftl_set(M,num)
  M.setoper('xftl',num.Value);
end

%-------------------------------------------------------------------------------
function oper_xftu_set(M,num)
  M.setoper('xftu',num.Value);
end

%-------------------------------------------------------------------------------
function oper_viscous_press(M,btn)
  M.setoper('visc',~M.oper.viscous);
  oper_viscous_set(M,btn);
end

%-------------------------------------------------------------------------------
function oper_viscous_set(M,btn)
  if (M.oper.viscous) btn.Text = 'Viscous'; btn.Value = true;
  else btn.Text = 'Inviscid'; btn.Value = false;
  end
end

%-------------------------------------------------------------------------------
function oper_givencl_press(M,btn)
  M.oper.givencl = ~M.oper.givencl;
  oper_givencl_set(M,btn);
end

%-------------------------------------------------------------------------------
function oper_givencl_set(M,btn)
  if (M.oper.givencl) btn.Text = 'Target cl'; btn.Value = true;
  else btn.Text = 'Target alpha'; btn.Value = false;
  end
end

%-------------------------------------------------------------------------------
function oper_initbl_press(M,btn)
  M.oper.initbl = ~M.oper.initbl;
  if (M.oper.initbl) btn.Text = 'Init BL'; btn.Value = false;
  else btn.Text = 'Reuse BL'; btn.Value = true;
  end
end

%-------------------------------------------------------------------------------
function oper_run_press(M,panels,axs,btn)
  pan_oper = panels{2};  gl_oper = pan_oper.Children(1); 
  pan_plot = panels{4};  gl_plot = pan_plot.Children(1); 
  M.solve();
  gl_oper.Children(2).Value = M.oper.alpha;
  M.oper.cltgt = M.post.cl;
  gl_oper.Children(4).Value = M.oper.cltgt;
  gl_plot.Children(2).Value = 1;
  gl_plot.Children(3).Data = get_post_table(M);
  set_conv_label(gl_oper.Children(17),M.glob.conv);
  plot_cp(M,axs);
  M.param.axplot = axs{2}; cla(axs{2}); axs{2}.NextPlot='replace'; 
  M.plot_airfoil(); M.plot_boundary_layer();
  V = axs{1}.InnerPosition; W = axs{2}.InnerPosition; 
  W(1)=V(1); W(3)=V(3); axs{2}.InnerPosition = W;
  plot_panels(M, panels);
end

%-------------------------------------------------------------------------------
function set_conv_label(lbl, conv)
  if (conv)
    lbl.Text = 'Converged'; lbl.FontColor = 'g'; lbl.FontWeight = 'normal';
    lbl.Tooltip = 'Coupled solver converged';
  else
    lbl.Text = 'NOT CONVERGED'; lbl.FontColor = 'r'; lbl.FontWeight = 'bold';
    lbl.Tooltip = 'Coupled solver did not converge; try running an easier case, e.g. a lower alpha, and reusing the BL solution.';
  end
end

%-------------------------------------------------------------------------------
function plot_cp(M,axs)
  M.param.axplot = axs{1}; cla(axs{1}); axs{1}.NextPlot='replace';  
  M.plot_cpplus();
end

%-------------------------------------------------------------------------------
function param_niglob_set(M,num)
  M.param.niglob = num.Value;
end

%-------------------------------------------------------------------------------
function param_rtol_set(M,num)
  M.param.rtol = num.Value;
end

%-------------------------------------------------------------------------------
function param_ncrit_set(M,num)
  M.param.ncrit = num.Value;
end

%-------------------------------------------------------------------------------
function param_verb_set(M,num)
  M.param.verb = num.Value;
end


%-------------------------------------------------------------------------------
function plot_dist_set(M,axs,dd)
  if (isempty(M.glob.U) && (dd.Value>1)), 
    f = msgbox('No viscous solution');
    dd.Value=1; 
    return; 
  end;
  switch dd.Value
    case 1
      plot_cp(M,axs);
    case 2
      plot_quantity(M, axs{1}, M.post.cf, 'c_f = skin friction coefficient');
    case 3
      plot_quantity(M, axs{1}, M.post.ds, '\delta^* = displacement thickness');
    case 4
      plot_quantity(M, axs{1}, M.post.th, '\theta = momentum thickness');
    case 5
      plot_quantity(M, axs{1}, M.post.Hk, 'H_k = kinematic shape parameter');
    case 6
      plot_quantity(M, axs{1}, M.post.ue, 'u_e = edge velocity');
    case 7
      plot_quantity(M, axs{1}, M.post.sa, 'amplification or c_{\tau}^{1/2}');
    case 8
      plot_quantity(M, axs{1}, M.post.Ret, 'Re_{\theta} = theta Reynolds number');
    otherwise
      disp('Unknown quantity for plotting');
  end
end

%-------------------------------------------------------------------------------
function plot_quantity(M, ax, q, qname)
% plots a quantity q over lower/upper/wake surfaces, on axes ax
  cla(ax); ax.NextPlot='replace'; 
  xy = [M.foil.x, M.wake.x];  % xy = M.isol.xi; % uncomment to plot vs xi
  ctype = {'r-', 'b-', 'k-'};
  if (isempty(M.vsol.Is)), plot(xy(1,:), q, [ctype{1},'o'], 'linewidth', 2);
  else
    sleg = {'lower', 'upper', 'wake'};
    for is = 1:3
      Is = M.vsol.Is{is};
      plot(ax,xy(1,Is), q(Is), [ctype{is},'o'], 'linewidth', 2, 'DisplayName', sleg{is});
      hold(ax,'on');
    end
  end
  set(ax, 'fontsize', 16); grid(ax, 'on');
  xlabel(ax, 'x = distance along chord', 'fontsize', 16);
  ylabel(ax, qname, 'fontsize', 16);
  if (~isempty(M.vsol.Is)), h = legend(ax,'Location', 'SouthEast'); set(h, 'fontsize', 14);  end
end

%-------------------------------------------------------------------------------
function T = get_post_table(M)
  names = {'lift coeff', 'moment coeff', 'drag coeff', 'skin fric drag coeff', 'press drag coeff'}';
  values = [M.post.cl, M.post.cm, M.post.cd, M.post.cdf, M.post.cdp]';
  svalues = cell(length(values),1);
  for i=1:length(values), svalues{i} = sprintf('%8.5f', values(i)); end
  T = table(names, svalues);
end


