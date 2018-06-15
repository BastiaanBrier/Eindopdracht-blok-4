Deze zip file bevat alle code nodig voor het runnen van het blast_script, insert_script en de webapplicatie.
Het blast script bevat een main functie die aangeroepen wordt als de code gerunt wordt en opent de files txt1 en txt2,
waar de geblaste forward en reverse sequenties in aanwezig zijn, om ze te blasten en ze met insert_script in de database
te stoppen. Aangezien alle meegeleverde sequenties uiteraard al in de database aanwezig zijn, zullen ze niet opnieuw in de
database gestopt worden.

De webapplicatie is lokaal te hosten door het webapp script te runnen, waarin de flask applicatie opgezet wordt. Hierin
wordt gebruik gemaakt van de HTML templates in de templates folder, die op hun beurt gebruik maken van de bootstrap toolkit
in de static folder.

De webapplicatie is naast de mogelijkheid hem lokaal te hosten ook benaderbaar via de volgende link:
http://cytosine.nl/~owe4_pg10/Project4/webapp.wsgi

Daar wordt via de wsgi file, waarvan tevens in deze zip een kopie te vinden is, naar webapp.py doorgelinkt.