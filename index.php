<!--https://www.it145.com/9/77315.html-->
<?php
mysql_connect("localhost","root","server@ntutrag");//與localhost連線、root是帳號、密碼處輸入自己設定的密碼
mysql_select_db("SRA_analysis");//我要從member這個資料庫撈資料
mysql_query("set names utf8");//設定utf8 中文字才不會出現亂碼
$data=mysql_query("select * from Final");//從member中選取全部(*)的資料
?>

<!DOCTYPE html>
<html lang="en">
<head>
<meta charset="UTF-8">

<!-- DataTables CSS -->
<link rel="stylesheet" type="text/css"
href="http://cdn.datatables.net/1.10.15/css/jquery.dataTables.css">
<!-- jQuery -->
<script type="text/javascript" charset="utf8" src="https://cdn.staticfile.org/jquery/1.10.2/jquery.min.js">
</script>
<!-- DataTables -->
<script type="text/javascript" charset="utf8"
src="http://cdn.datatables.net/1.10.15/js/jquery.dataTables.js"></script>


<title>SRA_analysis</title>
</head>

<body>
<div style="width:100%;">

<table id="mytable" border="1" style="width:100%;">
<thead>
<tr style="width:100%;">
<th style="width:10%;">Accession</th>
<th style="width:5%;">ST</th>
<th style="width:10%;">Serotype</th>
<th style="width:25%;">Inc type</th>
<th style="width:25%;">AMR</th>
<th style="width:25%;">Point</th>
</tr>
</thead>
<tbody>
<?php
for($i=1;$i<=mysql_num_rows($data);$i++)
{ $rs=mysql_fetch_row($data);
?>

<tr style="width:100%;">
<td style="width:10%; word-break:break-all;"><?php echo $rs[1]?></td>
<td style="width:5%; word-break:break-all;"><?php echo $rs[2]?></td>
<td style="width:15%; word-break:break-all;"><?php echo $rs[3]?></td>
<td style="width:20%; word-break:break-all;"><?php echo $rs[4]?></td>
<td style="width:25%; word-break:break-all;"><?php echo $rs[5]?></td>
<td style="width:25%; word-break:break-all;"><?php echo $rs[6]?></td>
</tr>

<?php }?>
</tbody>

</table>


<script>
$(document).ready( function () {
	$('#mytable').DataTable();
} );
</script>

</div>
</body>
</html>
