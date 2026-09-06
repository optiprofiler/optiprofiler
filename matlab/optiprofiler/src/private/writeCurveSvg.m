function writeCurveSvg(file, panels, heading)
%WRITECURVESVG Minimal, dependency-free visible output without MATLAB graphics.
%   This is deliberately a plain line chart, not a replacement graphics engine.
%   Curves are 2-by-N numeric arrays; nonfinite points break a line rather than
%   inventing data. Callers apply documented display-only history protection.
    columns = min(2, max(1, numel(panels)));
    rows = max(1, ceil(numel(panels) / columns));
    legend_rows = 0;
    for k=1:numel(panels)
        count=0;
        for j=1:numel(panels(k).labels), count=count+numel(wrapLegend(panels(k).labels{j})); end
        legend_rows=max(legend_rows,count);
    end
    panel_height = 430 + 16*legend_rows;
    width = 680 * columns; height = panel_height * rows + 65;
    fid = fopen(file, 'w');
    if fid < 0, error('OptiProfiler:SvgOutput', 'Cannot write ''%s''.', file); end
    cleanup = onCleanup(@() closeIfOpen(fid));
    fprintf(fid, '<svg xmlns="http://www.w3.org/2000/svg" width="%d" height="%d" viewBox="0 0 %d %d">\n', width, height, width, height);
    fprintf(fid, '<rect width="100%%" height="100%%" fill="white"/>\n');
    textAt(fid, 24, 26, heading, 18, 'start');
    textAt(fid, 24, 49, 'Portable SVG fallback; numeric results and profile scores are unchanged.', 12, 'start');
    colors = {'#0072b2', '#d55e00', '#009e73', '#cc79a7', '#e69f00', '#56b4e9', '#333333'};
    for k = 1:numel(panels)
        p = panels(k);
        ox = mod(k-1, columns)*680; oy = floor((k-1)/columns)*panel_height+65;
        left = ox+90; top = oy+70; w = 545; h = 220;
        textAt(fid, ox+24, oy+24, p.title, 15, 'start');
        xs = []; ys = [];
        for j = 1:numel(p.curves)
            c = p.curves{j};
            if isempty(c), continue; end
            valid = all(isfinite(c), 1);
            xs = [xs, c(1,valid)]; ys = [ys, c(2,valid)]; %#ok<AGROW>
        end
        if isfield(p,'include_zero') && p.include_zero, ys=[ys,0]; end
        [xmin,xmax,xscale] = limits(xs); [ymin,ymax,yscale] = limits(ys);
        fprintf(fid, '<rect x="%.5g" y="%.5g" width="%g" height="%g" fill="none" stroke="#888"/>\n', left, top, w, h);
        if isfield(p,'include_zero') && p.include_zero
            zero_y=top+h*(1-(-ymin/yscale)/(ymax/yscale-ymin/yscale));
            fprintf(fid,'<line x1="%g" x2="%g" y1="%.5g" y2="%.5g" stroke="#aaa" stroke-dasharray="4 3"/>\n',left,left+w,zero_y,zero_y);
        end
        for t = 0:4
            fraction = t/4;
            textAt(fid, left+w*fraction, top+h+20, sprintf('%.4g', (1-fraction)*xmin+fraction*xmax), 11, 'middle');
            textAt(fid, left-8, top+h*(1-fraction)+4, sprintf('%.4g', (1-fraction)*ymin+fraction*ymax), 11, 'end');
        end
        textAt(fid, left+w/2, top+h+43, p.xlabel, 12, 'middle');
        textAt(fid, left, top-12, p.ylabel, 12, 'start');
        legend_row=0;
        for j = 1:numel(p.curves)
            c = p.curves{j};
            color = colors{mod(j-1,numel(colors))+1};
            if ~isempty(c)
                valid = all(isfinite(c), 1);
                boundaries = diff([false,valid,false]);
                starts = find(boundaries==1); ends = find(boundaries==-1)-1;
                for q = 1:numel(starts)
                    index = starts(q):ends(q);
                    % Normalize before subtracting so even opposite-signed
                    % finite values near realmax cannot overflow the range.
                    x = left+w*(c(1,index)/xscale-xmin/xscale)/(xmax/xscale-xmin/xscale);
                    y = top+h*(1-(c(2,index)/yscale-ymin/yscale)/(ymax/yscale-ymin/yscale));
                    fprintf(fid, '<polyline fill="none" stroke="%s" stroke-width="1.7" points="', color);
                    fprintf(fid, '%.5g,%.5g ', [x;y]); fprintf(fid, '"/>\n');
                    if numel(index)==1
                        fprintf(fid,'<circle cx="%.5g" cy="%.5g" r="3" fill="%s"/>\n',x,y,color);
                    end
                end
            end
            % Every legend entry is plain escaped text; no LaTeX installation
            % or interpretation is assumed by this fallback.
            legend_x = ox+24; legend_y = oy+405+legend_row*16;
            fprintf(fid, '<rect x="%g" y="%g" width="12" height="3" fill="%s"/>', legend_x, legend_y-5, color);
            lines=wrapLegend(p.labels{j});
            for line=1:numel(lines)
                textAt(fid,legend_x+17,legend_y+(line-1)*16,lines{line},10,'start');
            end
            legend_row=legend_row+numel(lines);
        end
        if isempty(xs), textAt(fid,left+w/2,top+h/2,'No finite curve points',13,'middle'); end
        if isfield(p,'note') && ~isempty(p.note)
            notes=splitlines(string(p.note));
            for j=1:numel(notes), textAt(fid,ox+24,oy+351+12*j,notes(j),10,'start'); end
        end
    end
    fprintf(fid, '</svg>\n');
    message=ferror(fid); closed=fclose(fid); clear cleanup;
    if ~isempty(message) || closed~=0
        error('OptiProfiler:SvgOutput','Could not finish writing ''%s'': %s',file,message);
    end
end

function lines=wrapLegend(value)
    % Fixed plain-text wrapping also works without a font measurement engine.
    % Sixty characters at 10 px fit even wide Latin glyphs in this panel.
    value=char(string(value)); value=strrep(value,newline,' ');
    starts=1:60:max(1,numel(value)); lines=cell(1,numel(starts));
    for k=1:numel(starts), lines{k}=value(starts(k):min(end,starts(k)+59)); end
end

function [lo,hi,scale] = limits(values)
    if isempty(values), lo=0; hi=1; scale=1; return; end
    lo=min(values); hi=max(values);
    if lo==hi
        if lo>0, lo=0.95*lo;
        elseif hi<0, hi=0.95*hi;
        else, lo=-1; hi=1;
        end
        if lo==hi
            % A five-percent margin can round back to the same subnormal.
            if lo>0, lo=0; else, hi=0; end
        end
    end
    scale=max(abs([lo,hi]));
end

function textAt(fid,x,y,text,size,anchor)
    text=char(string(text));
    text(double(text)<32 & ~ismember(double(text),[9,10,13]))='?';
    text=strrep(strrep(strrep(text,'&','&amp;'),'<','&lt;'),'>','&gt;');
    text=strrep(strrep(text,'"','&quot;'),'''','&apos;');
    fprintf(fid,'<text x="%g" y="%g" font-family="sans-serif" font-size="%g" text-anchor="%s">%s</text>\n',x,y,size,anchor,text);
end

function closeIfOpen(fid)
    if ~isempty(fopen(fid)), fclose(fid); end
end
