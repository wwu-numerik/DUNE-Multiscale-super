| Image  | Info |
| :----- | :--- |
define(`slug', `BASENAME-$1')dnl
define(`SVG', `https://images.microbadger.com/badges/image/dunecommunity/slug($1).svg')dnl
define(`link', `@<:@!@<:@@>:@(SVG($1))@>:@')dnl
define(`section',
`| slug($1) | link($1)(http://microbadger.com/images/dunecommunity/slug($1) "slug($1) layer") |
')dnl
include(m4sugar.m4)
# output to stdout
m4_divert_push(1)dnl
m4_foreach_w(REPO, REPONAMES, section(REPO))
