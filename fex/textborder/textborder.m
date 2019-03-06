function textborder(ax, x, y, string, text_color, border_color, varargin)
%TEXTBORDER Display text with border.
%   TEXTBORDER(AX, X, Y, STRING)
%   Creates text on the current figure with a one-pixel border around it.
%   The default colors are white text on a black border, which provides
%   high contrast in most situations.
%   
%   TEXTBORDER(AX, X, Y, STRING, TEXT_COLOR, BORDER_COLOR)
%   Optional TEXT_COLOR and BORDER_COLOR specify the colors to be used.
%   
%   Optional properties for the native TEXT function (such as 'FontSize')
%   can be supplied after all the other parameters.
%   Since usually the units of the parent axes are not pixels, resizing it
%   may subtly change the border of the text out of position. Either set
%   the right size for the figure before calling TEXTBORDER, or always
%   redraw the figure after resizing it.
%   
%   Author: Jo�o F. Henriques, April 2010
%   Edit by Geoff Stanley, 09/05/2018, to pass ax as first argument

	if isempty(string), return; end
	
	if nargin < 6, border_color = 'k'; end  %default: black border
	if nargin < 5, text_color = 'w'; end  %default: white text
	
	%border around the text, composed of 4 text objects
	offsets = [0 -1; -1 0; 0 1; 1 0];
	for k = 1:4
		h = text(ax, x, y, string, 'Color',border_color, varargin{:});
		
		%add offset in pixels
		set(h, 'Units','pixels')
		pos = get(h, 'Position');
		set(h, 'Position', [pos(1:2) + offsets(k,:), 0])
		set(h, 'Units','data')
	end
	
	%the actual text inside the border
	h = text(ax, x, y, string, 'Color',text_color, varargin{:});
	
	%same process as above but with 0 offset; corrects small roundoff
	%errors
	set(h, 'Units','pixels')
	pos = get(h, 'Position');
	set(h, 'Position', [pos(1:2), 0])
	set(h, 'Units','data')

end

