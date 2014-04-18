/* 
This file is part of the COPS web page.

 COPS web page by Miguel Ferreira, Nuno Roma, and Luis M. S. Russo / KDBIO /
INESC-ID is licensed under a Creative Commons Attribution 3.0 Unported
License.  Permissions beyond the scope of this license may be available at #email.

 Software is distributed in the hope that it will be useful, but
WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.

 License available at http://creativecommons.org/licenses/by/3.0/
*/

var sections=["email","down","doc","home"]; // Array of sections

function setSections()
{
//    alert(location.hash);
    for (x in sections)
    {
        if("#"+sections[x] == location.hash)
        {
            document.getElementById(sections[x]).className="activeArea";
            document.getElementById(sections[x]+"Txt").style.visibility="visible";
        }
        else
        {
            document.getElementById(sections[x]).className="area";
            document.getElementById(sections[x]+"Txt").style.visibility="hidden";
        }
    }
}

function setHash(nh)
{
    window.location.assign("#"+nh)
    setSections();
}

function load()
{
    if(location.hash=="")
        setHash('home');
    setSections();
/* Set the onclick functions */    
    for (x in sections)
        document.getElementById(sections[x]).onclick=Function("setHash('"+sections[x]+"')");
}




