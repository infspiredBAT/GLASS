jsDeleteRow <- "shinyjs.delRow = function(id)
{$('.samplesdt tr').eq(id).hide(); }"

jsSwapRow <- "shinyjs.swapRow = function(id)
{
console.log('swap'+id);
var swp = $('.samplesdt tr:eq('+ id +') td:eq(1)').text();
console.log(swp);
$('.samplesdt tr:eq('+ id +') td:eq(1)').text($('.samplesdt tr:eq('+ id +') td:eq(3)').text());
$('.samplesdt tr:eq('+ id +') td:eq(3)').text(swp)}"

