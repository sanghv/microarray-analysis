I use MikTex on Windoze so things are pretty straight forward to getting presentation.tex to compile (automatically download packages that are missing).


I made sure to use a default beamer theme so it should be as simple as downloading the beamer package automatically:

\documentclass[serif]{beamer}

(IF IT DOESN'T DOWNLOAD AUTOMATICALLY, I HATE TO SAY IT BUT JUST SCOUR THE INTERNET FOR INSTRUCTIONS OR GO RIGHT TO THEIR BITBUCKET PAGE https://bitbucket.org/rivanvx/beamer/wiki/Home)

These packages should be able to be downloaded automatically:

\usepackage{amsmath}
\usepackage{graphicx}
\usepackage{rotating}
\usepackage{etoolbox}
\usepackage{xfrac} % for \sfrac{...}{...}
\usepackage{array} % for centered tabular
\usepackage{color}
\usepackage{listings}

The layout of the beamer project is done so the only thing left is frames


use \begin{frame}  <CONTENT HERE \end{frame} to make a new frame

TO MAKE A FRAME TITLE:
\frametitle{ <TITLE HERE> }
TO MAKE A FRAME SUBTITLE:
\framesubtitle{ <SUBTITLE HERE> }

TO MAKE A BLOCK:
\begin{block}{ <BLOCK TITLE> }
  <CONTENT HERE>
\end{block}


TO MAKE A BLOCK WITHOUT A TITLE:
\begin{block}{ \vspace{-0.5in} }
	<CONTENT HERE>
\end{block}



I will take care of the bibliography. You can just make notations using a standard Latex comment on whether you cited anything:

i.e.
% CITE: <INSERT SOME LINK>


HOW DO I MAKE ALL THOSE TRANSITION THINGS?!
Please see this guide or the plentiful examples in our presentation to understand how to do split animations over multiple slides:
http://www.math-linux.com/spip.php?article77   (AWESOME TUTORIAL)


Don't worry too much about the code formatting. I used the listing package to do some of the stuff. I don't know a lot about it and just do it as I go. Just google and use stackoverflow if you can't get it to do something.
