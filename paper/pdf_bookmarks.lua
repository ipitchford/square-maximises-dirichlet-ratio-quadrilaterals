-- PDF outline only: math in section titles gets a plain Unicode alternative (\texorpdfstring); the typeset page is unchanged.
local subs = {
  {"\\mathrm{([^}]*)}", "%1"}, {"\\theta", "θ"}, {"\\to", "→"}, {"\\le", "≤"}, {"\\pi", "π"},
  {"%^\\circ", "°"}, {"\\circ", "°"}, {"'", "′"}, {"_", ""}, {"%$", ""}, {"\\,", " "},
}
local function plain(tex)
  local s = tex
  for _, p in ipairs(subs) do s = s:gsub(p[1], p[2]) end
  s = s:gsub("\\", "")
  return s
end
function Header(h)
  h.content = h.content:walk({
    Math = function(m)
      -- typeset: the math itself (raw, so pandoc leaves it out of its bookmark string);
      -- bookmark: a plain string that \texorpdfstring{}{...} hides from the typeset page
      return {pandoc.RawInline("latex", "\\(" .. m.text .. "\\)"), pandoc.RawInline("latex", "\\texorpdfstring{}{"),
              pandoc.Str(plain(m.text)), pandoc.RawInline("latex", "}")}
    end
  })
  return h
end
