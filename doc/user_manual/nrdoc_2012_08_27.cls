%%
%% This is file nrdoc.cls
%%
%% Based on nrdoc.cls by
%% Harald H. Soleng <harald.soleng@nr.no>
%% Anders Løland <anders.loland@nr.no>
%%
%% by Lars Brusletto Sveen <lars@munch-design.no>
%%
%% Edited by H. H. Soleng, February-May 2005

\NeedsTeXFormat{LaTeX2e}[1995/12/01]

\ProvidesClass{nrdoc}[2007/10/08 v2.2.11-1 NR report class]

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Adding logic and math to latex for easy programming
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\RequirePackage{ifthen}             % for logic in class file
\RequirePackage{calc}               % for math in class file
\RequirePackage{textcase}
\RequirePackage{ifpdf}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Adding new variables
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\ifpdf
   \pdfcompresslevel=9
\fi

%\newboolean{pdfdataset}
%\setboolean{pdfdataset}{false}

\newboolean{aboutauthors}\setboolean{aboutauthors}{false}
\newboolean{qualityassurance}\setboolean{qualityassurance}{false}

\newboolean{longreport}
\newboolean{onesided} 
\newboolean{draftversion}
\newboolean{report}
\newboolean{screen}
\newboolean{showlabels}
\newboolean{closedavailability}
\newboolean{noindentpar}
\setboolean{noindentpar}{true}
\newboolean{indentpar}
\setboolean{indentpar}{false}
\newboolean{firstindent}
\newboolean{language_selected}
\newboolean{norsk_selected}
\newboolean{nynorsk_selected}
\newboolean{english_selected}
\newboolean{citenumeric}
\newboolean{citealphanumeric}
\newboolean{hidden}
\newboolean{dvips}
\newboolean{nodesign}
\newboolean{noindex}
\newboolean{nobibstyle}
\newboolean{utf8}
\setboolean{utf8}{false}

\setboolean{nobibstyle}{false}
\setboolean{noindex}{false}


\newcounter{numberofauthors}
\setcounter{numberofauthors}{1}

\newcounter{largenumberofauthors}
\setcounter{largenumberofauthors}{4}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Adding new functions
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\providecommand{\authorname}            {Author\ifthenelse{\value{numberofauthors}>1}{s}{}}
\providecommand{\availabilityname}                      {Availability}
\providecommand{\availabilitylevel}                     {Confidential}
\providecommand{\availabilitylevelopen}                 {Open}
\providecommand{\availabilitylevelinternal}             {Internal}
\providecommand{\availabilitylevelconfidential}         {Confidential}
\providecommand{\availabilitylevelstatoilhydro}         {StatoilHydro restricted}
\providecommand{\availabilitylevelstatoil}         {Statoil restricted}
\providecommand{\availabilitylevelstrictlyconfidential} {Strictly confidential}
\providecommand{\closedavailable}       {Closed}
\providecommand{\datename}              {Date}
\providecommand{\qaname}              	{Quality assurance}
\providecommand{\examplename}           {Example}
\providecommand{\exercisename}          {Exercise}
\providecommand{\emailname}             {e-mail}
\providecommand{\keywordsname}          {Keywords}
\providecommand{\listexamplename}       {List of Examples}
\providecommand{\listexercisename}    	{List of Exercises}
\providecommand{\numberofpagesname}     {Number of pages}
\providecommand{\notename}              {Note}
\providecommand{\nrdocument}            {NR Document}
\providecommand{\openavailable}         {Open}
\providecommand{\projectname}           {Project}
\providecommand{\projectnumbername}     {Project number}
\providecommand{\publicationnumbername} {Publication number}
\providecommand{\reportnumbername}      {Report no}
\providecommand{\notenumbername}        {Note no}
\providecommand{\researchfieldname}     {Research field}
\providecommand{\reportname}            {Report}
\providecommand{\targetgroupname}       {Target group}
\providecommand{\titlename}             {Title}
\providecommand{\yearname}              {Year}
\providecommand{\isbnname}              {ISBN}
\providecommand{\nrname}                {Norwegian Computing Center}
\providecommand{\aboutauthorsheading}   {The author\ifthenelse{\value{numberofauthors}>1}{s}{}}


\providecommand{\aboutauthors}[1]{\setboolean{aboutauthors}{true}\renewcommand{\aboutauthorstext}{#1}}
\providecommand{\qa}[1]{\setboolean{qualityassurance}{true}\renewcommand{\printqa}{#1}}

\providecommand{\designelement}[1]{%
  \ifthenelse{\boolean{nodesign}}{%
    \includegraphics[draft]{#1}
  }{%
    \includegraphics{#1}        
  }
}

\providecommand{\email}[1]{%
  \ifthenelse{\boolean{titlepage}}{%
    \ifthenelse{\equal{#1}{}}{}{%
      \xspace\href{mailto:#1}{\texttt{<#1>}}}%
  }%
  {\ignorespaces}%
  \renewcommand{\and}{%
     \ifthenelse{\boolean{hidden}}%
        {\addtocounter{numberofauthors}{1}}
  }     
  \newline%     
}





\providecommand{\printyear}{\number\year}
\providecommand{\publicationyear}[1]{\renewcommand{\printyear}{#1}}

\providecommand{\printshorttitle}{\@title}
\providecommand{\shorttitle}[1]{\renewcommand{\printshorttitle}{#1}}

\providecommand{\@subtitle}{}
\providecommand{\printsubtitle}{\@subtitle}
\providecommand{\subtitle}[1]{\renewcommand{\@subtitle}{#1}}

\providecommand{\printtitle}{\@title}
\providecommand{\nrtitle}[1]{\renewcommand{\printtitle}{#1}}

\providecommand{\printisbn}{\typeout{NRdoc warning: ISBN number not set}}
\providecommand{\isbn}[1]{%
  \renewcommand{\printisbn}{#1}}
\providecommand{\shortisbn}[1]{%
  \renewcommand{\printisbn}{82-536-#1}}

\providecommand{\printreportnumber}{%
  \typeout{NRdoc warning: Report number not set}No number}
\providecommand{\reportnumber}[1]{\renewcommand{\printreportnumber}{%
    \MakeUppercase{#1}}}
\providecommand{\pdfprintkeywords}{%
}
\providecommand{\printqa}{%
  \typeout{NRdoc warning: QA person undefined}Unknown}
\providecommand{\printqatext}{%
  \ifthenelse{\boolean{qualityassurance}}{%
       \titleitem{\qaname}	\printqa}{}}
\providecommand{\printkeywords}{%
  \protect\typeout{NRdoc warning: Keywords not defined}}
\providecommand{\keywords}[1]{\renewcommand{\printkeywords}{#1}%
                              \renewcommand{\pdfprintkeywords}{#1}}

\providecommand{\printsubject}{}
\providecommand{\subject}[1]{\renewcommand{\printsubject}{#1}}

\providecommand{\printauthorforpdf}{Norsk Regnesentral}
\providecommand{\authorpdf}[1]{\renewcommand{\printauthorforpdf}{#1}}

\providecommand{\printproject}{%
                \typeout{NRdoc warning: Project not defined}}
\providecommand{\project}[1]{\renewcommand{\printproject}{#1}}

\providecommand{\printprojectnumber}{%
\typeout{NRdoc warning: Project number not defined}}
\providecommand{\projectnumber}[1]{\renewcommand{\printprojectnumber}{#1}}

\providecommand{\printtarget}{%
                \typeout{NRdoc warning: Target group not defined}}
\providecommand{\target}[1]{\renewcommand{\printtarget}{#1}}

\providecommand{\printavailability}{%
                \typeout{NRdoc warning: Availability level not
defined}}

\providecommand{\availability}[1]{%
\ifthenelse{\equal{#1}{Open}}{\renewcommand{\availabilitylevel}{\availabilitylevelopen}}{%
\ifthenelse{\equal{#1}{Internal}}{\renewcommand{\availabilitylevel}{\availabilitylevelinternal}}{%
\ifthenelse{\equal{#1}{StatoilHydro restricted}}{\renewcommand{\availabilitylevel}{\availabilitylevelstatoilhydro}}{%
\ifthenelse{\equal{#1}{Statoil restricted}}{\renewcommand{\availabilitylevel}{\availabilitylevelstatoil}}{%
\ifthenelse{\equal{#1}{Confidential}}{}{%
\ifthenelse{\equal{#1}{Strictlyconfidential}}{\renewcommand{\availabilitylevel}{\availabilitylevelstrictlyconfidential}}{%
\typeout{NRdoc warning: Unknown availability level: #1} %
}}}}}}
\renewcommand{\printavailability}{\availabilitylevel}}

\providecommand{\printresearchfield}{%
\typeout{NRdoc warning: Reserach field not defined}}
\providecommand{\researchfield}[1]{\renewcommand{\printresearchfield}{#1}}

\newboolean{frontpagefiguredefined}
\setboolean{frontpagefiguredefined}{false}

\providecommand{\frontpagefigurefile}{}
\providecommand{\frontpagefigurecpright}{}
\providecommand{\frontpagefigure}[2][{}]{%
\renewcommand{\frontpagefigurefile}{#2}%
\renewcommand{\frontpagefigurecpright}{#1}%
\setboolean{frontpagefiguredefined}{true}}

\providecommand{\frontpagefigurecommands}{}
\providecommand{\frontpagefigurecommand}[1]{%
\renewcommand{\frontpagefigurecommands}{#1}}

\renewcommand{\and}{\ifthenelse{\boolean{hidden}}{\addtocounter{numberofauthors}{1}}{\linebreak}}

\newboolean{titlepage}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Handling options
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\DeclareOption{noindex}{\setboolean{noindex}{true}}
\DeclareOption{utf8}{\setboolean{utf8}{true}}

\DeclareOption{long}{\setboolean{longreport}{true}}
\DeclareOption{10pt}{%
\PassOptionsToClass{10pt}{book}\PassOptionsToClass{10pt}{article}}
\DeclareOption{11pt}{%
\PassOptionsToClass{11pt}{book}\PassOptionsToClass{11pt}{article}}
\DeclareOption{12pt}{%
\PassOptionsToClass{12pt}{book}\PassOptionsToClass{12pt}{article}}

\DeclareOption{draft}{%
\setboolean{draftversion}{true}%
\setboolean{showlabels}{true}%
\PassOptionsToPackage{draft}{graphicx}%
\PassOptionsToClass{draft}{book}%
\PassOptionsToClass{draft}{article}%
}

\DeclareOption{showlabels}{\setboolean{showlabels}{true}}
\DeclareOption{nobibstyle}{\setboolean{nobibstyle}{true}}

\DeclareOption{note}{\setboolean{report}{false}}
\DeclareOption{report}{\setboolean{report}{true}}


\DeclareOption{dvips}{\setboolean{dvips}{true}}

\DeclareOption{nodesign}{\setboolean{nodesign}{true}}

\DeclareOption{screen}{\setboolean{screen}{true}%
  \setboolean{onesided}{true}%    
  \typeout{NRdoc: Setting onesided true}%
  \PassOptionsToClass{oneside}{book}%
  \PassOptionsToClass{oneside}{article}
}

\DeclareOption{oneside}{%
  \setboolean{onesided}{true}%
  \typeout{NRdoc: Setting onesided true}%
  \PassOptionsToClass{oneside}{book}%
  \PassOptionsToClass{oneside}{article}
}%

\DeclareOption{twoside}{% 
  \setboolean{onesided}{false}% 
  \PassOptionsToClass{twoside}{book}%
  \PassOptionsToClass{twoside}{article}%
  \ifthenelse{\boolean{screen}}{%
  \ClassError{nrdoc}{%
The options `twoside' and `screen' are incompatible}{%
Oh dear! Have you ever tried twosided printing on a %
computer screen?\MessageBreak 
It could damage your eyes. I rather advise you to
remove either\MessageBreak the `screen' or the `twoside' option in %
the documentclass command!\MessageBreak}}{}
}%



\DeclareOption{closed}{\setboolean{closedavailability}{true}}
\DeclareOption{noindentpar}{\setboolean{noindentpar}{true}}
\DeclareOption{indentpar}{\setboolean{noindentpar}{false}}
% The bluelines option is obsolete (kept for backwards compatability)
\DeclareOption{bluelines}{\typeout{NRdoc warning: Option bluelines is obsolete.}}
\DeclareOption{american}{\PassOptionsToPackage{american}{babel}
                         \setboolean{language_selected}{true}   
                         \setboolean{firstindent}{false}
                         \setboolean{english_selected}{true}}
\DeclareOption{british}{\PassOptionsToPackage{british}{babel}
                        \setboolean{language_selected}{true}    
                        \setboolean{firstindent}{false}
                        \setboolean{english_selected}{true}}
\DeclareOption{english}{\PassOptionsToPackage{english}{babel}
                        \setboolean{language_selected}{true}    
                         \setboolean{firstindent}{false}
                        \setboolean{english_selected}{true}}
\DeclareOption{UKenglish}{\PassOptionsToPackage{UKenglish}{babel}
                          \setboolean{language_selected}{true}
                          \setboolean{firstindent}{false}
                          \setboolean{english_selected}{true}}
\DeclareOption{USenglish}{\PassOptionsToPackage{USenglish}{babel}
                          \setboolean{language_selected}{true}  
                          \setboolean{firstindent}{false}
                          \setboolean{english_selected}{true}}
\DeclareOption{norsk}{\PassOptionsToPackage{norsk}{babel}
                      \setboolean{norsk_selected}{true} 
                      \setboolean{firstindent}{true}
                      \setboolean{language_selected}{true}}
\DeclareOption{nynorsk}{\PassOptionsToPackage{nynorsk}{babel}
                        \setboolean{nynorsk_selected}{true}     
                        \setboolean{firstindent}{true}
                        \setboolean{language_selected}{true}}
% The times option is obsolete (kept for backwards compatability)
\DeclareOption{times}{\typeout{NRdoc warning: Option times is obsolete.}}




\DeclareOption{citenumeric}{%
  \setboolean{citenumeric}{true}}
\DeclareOption{citealphanumeric}{%
  \setboolean{citenumeric}{true}
  \setboolean{citealphanumeric}{true}
}

\DeclareOption*{\PackageWarning{nrdoc}{Option `\CurrentOption' ignored.}}

\ExecuteOptions{10pt,note,twoside} 


\ProcessOptions\relax


%  Norwegian translation of cover and title side expressions
\ifthenelse{\boolean{norsk_selected} \or
            \boolean{nynorsk_selected}}{%
        \renewcommand{\authorname}              {Forfatter\ifthenelse{\value{numberofauthors}>1}{e}{}}
        \renewcommand{\availabilityname}        {Tilgjengelighet}
        \renewcommand{\availabilitylevel}                       {Konfidensiell}
        \renewcommand{\availabilitylevelopen}                   {{\AA}pen}
        \renewcommand{\availabilitylevelinternal}               {Intern}
        \renewcommand{\availabilitylevelconfidential}           {Konfidensiell}
        \renewcommand{\availabilitylevelstatoilhydro}           {StatoilHydro begrenset}
        \renewcommand{\availabilitylevelstatoil}           {Statoil begrenset}
        \renewcommand{\availabilitylevelstrictlyconfidential}   {Strengt konfidensiell}
        \renewcommand{\datename}                {Dato}
	\renewcommand{\qaname}              	{Kvalitetssikring} 
        \renewcommand{\examplename}             {Eksempel}
        \renewcommand{\exercisename}            {Oppgave}
        \renewcommand{\emailname}               {e-post}
        \renewcommand{\keywordsname}            {Emneord}
        \renewcommand{\listexamplename}         {Eksempler}
        \renewcommand{\listexercisename}        {Oppgaver}
        \renewcommand{\numberofpagesname}       {Antall sider}
        \renewcommand{\notename}                {Notat}
        \renewcommand{\openavailable}           {{\AA}pen}
        \renewcommand{\reportname}              {Rapport}
        \renewcommand{\projectname}             {Prosjekt}
        \renewcommand{\projectnumbername}       {Prosjektnummer}
        \renewcommand{\publicationnumbername}   {Publikasjonsnummer}
        \renewcommand{\reportnumbername}        {Rapportnr}
        \renewcommand{\notenumbername}          {Notatnr}
        \renewcommand{\researchfieldname}       {Satsningsomr\aa de}
        \renewcommand{\targetgroupname}         {M\aa lgruppe}
        \renewcommand{\titlename}               {Tittel}
        \renewcommand{\yearname}                {{\AA}r}
        \ifthenelse{\boolean{longreport}}{
        \providecommand{\abstractname}          {Sammendrag}}{}
        \renewcommand{\nrname}                  {Norsk Regnesentral}
        \renewcommand{\aboutauthorsheading}     {Forfatter\ifthenelse{\value{numberofauthors}>1}{ne}{en}}

%  Nynorsk
\ifthenelse{\boolean{nynorsk_selected}}{%
        \renewcommand{\listexamplename}         {Eksempel}
        \renewcommand{\listexercisename}        {Oppg\aa ver}
        \renewcommand{\exercisename}        	{Oppg\aa ve}
        \renewcommand{\authorname}              {Forfattar\ifthenelse{\value{numberofauthors}>1}{ar}{}}
        \renewcommand{\availabilityname}        {Tilgjengeleg}
        \renewcommand{\availabilitylevelopen}   {Open}
        \renewcommand{\numberofpagesname}       {Tal p\aa\ sider}
        \ifthenelse{\boolean{longreport}}{
        \providecommand{\abstractname}          {Samandrag}}{}
        \renewcommand{\aboutauthorsheading}     {Forfattar\ifthenelse{\value{numberofauthors}>1}{ane}{en}}
        }{}
}{\ifthenelse{\boolean{longreport}}{
        \providecommand{\abstractname}          {Abstract}}{}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Loading other packages and the base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\RequirePackage[T1]{fontenc}

\RequirePackage{amsmath}        % for advanced math
\RequirePackage{fancyhdr}       % for headers
\RequirePackage{rotating}       % for text rotation
\RequirePackage{makeidx}        % for index
\RequirePackage{multicol}       % for multiple columns
\RequirePackage{paralist}       % for list environments
\RequirePackage{moreverb}       % extra verbatims
\RequirePackage{xspace}         % for space control
\RequirePackage{textcomp}       % for special characters
\RequirePackage{array}          % for tabular functionality
\RequirePackage{float}          % for float control
\RequirePackage[plain]{algorithm}       % algorithm



\ifpdf\RequirePackage{thumbpdf} \fi      % for thumbnails

\RequirePackage{afterpage}      % for float control
\RequirePackage{color}          % for colour
\RequirePackage{graphicx}
      
\ifthenelse{\boolean{language_selected}}{%
  \RequirePackage{babel}
}{}

\ifthenelse{\boolean{utf8}}{%
\RequirePackage[utf8]{inputenc}}{%
\RequirePackage[latin1]{inputenc}}



  \ifthenelse{\boolean{longreport}}%
  {\providecommand{\theHalgorithm}{\thechapter.\arabic{algorithm}}%
   \LoadClass{book}[2000/05/19]}%
  {\providecommand{\theHalgorithm}{\arabic{algorithm}}%
   \LoadClass{article}[2000/05/19]}

\def\toclevel@algorithm{0}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Loading other packages after the base class
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


\ifthenelse{\boolean{showlabels}}{%
 \RequirePackage{showlabels}}

\ifthenelse{\boolean{firstindent}}{% 
 \RequirePackage{indentfirst}}{}

% Set Palatino as the text font. Mathpazo loads Palatino-like fonts
% for math mode.
%
% For now, set Helvetica as the sans serif fonts. This should be
% changed to Arial in a future version.
\renewcommand{\rmdefault}{ppl}
\renewcommand{\sfdefault}{phv}
\renewcommand{\ttdefault}{aett}

% Removes the much hated Palatino math font
%\RequirePackage{mathpazo}

% The leading should be increased slightly for Palatino. This setting
% gives a leading of approx. 14pt for 10 pt text.
\renewcommand{\baselinestretch}{1.167}

\ifthenelse{\boolean{citenumeric}}{%  
 \RequirePackage[numbers,sort&compress]{natbib}   
}{%
 \RequirePackage[sort,round]{natbib}
}


% NR profile blue in CMYK values
\definecolor{nrblue}{cmyk}{1,0.66,0,0}

\ifthenelse{\boolean{screen}}{
% use hypertex tools    
\@twosidefalse
\ifpdf
  \RequirePackage[pdftex,colorlinks=true,pdfpagelabels,%
                  plainpages=false,citecolor=nrblue,%
                  backref=page]{hyperref}
\else
  \ifthenelse{\boolean{dvips}}{
  \RequirePackage[dvips,colorlinks=true,pdfpagelabels,
                  plainpages=false,citecolor=nrblue,%
                  backref=page]{hyperref}}{%
  \RequirePackage[dvipdfm,colorlinks=true,pdfpagelabels,
                  plainpages=false,citecolor=nrblue,%
                  backref=page]{hyperref}
  }
\fi
  }{%
% Set all hypertex links black %
\ifpdf
  \RequirePackage[pdftex,colorlinks=true,pdfpagelabels,plainpages=false,%
                  linkcolor=black,citecolor=black,%
                  pagecolor=black,urlcolor=black,bookmarks]{hyperref}
\else 
  \ifthenelse{\boolean{dvips}}{
  \RequirePackage[dvips,colorlinks=true,pdfpagelabels,plainpages=false,%
                  linkcolor=black,citecolor=black,%
                  pagecolor=black,urlcolor=black]{hyperref}}{%
  \RequirePackage[dvipdfm,colorlinks=true,pdfpagelabels,plainpages=false,%
                  linkcolor=black,citecolor=black,%
                  pagecolor=black,urlcolor=black]{hyperref}
  }
\fi
}



% Load parskip to get ``block'' paragraphs and set parskip. Note that
% this may have adverse effects since \parskip is used in many places.
% The package corrects the behaviour in the toc and in lists. Spacing
% around sections is corrected through titlesec. Further corrections
% may be needed.
\ifthenelse{\boolean{noindentpar}}{%
  \typeout{No indents in paragraphs}
  \RequirePackage{parskip}[2001/04/09]
  \parskip=0.5\baselineskip
  \advance\parskip by 0pt plus 2pt
}{%
  \typeout{Using paragraph indentation}
    
}

% Load geometry and set page layout
\RequirePackage{geometry}[2002/07/08]

\ifthenelse{\boolean{screen}}{%
\geometry{a4paper, left=30mm, right=30mm, 
top=20mm, bottom=35mm,
footskip=20mm, footnotesep=\baselineskip, nohead, twoside=false}}{}

\ifthenelse{\boolean{onesided}}{%
\typeout{Selecting onesided page geometry}
\geometry{a4paper, left=30mm, right=30mm, top=20mm, bottom=35mm,
  footskip=20mm, footnotesep=\baselineskip, nohead, twoside=false}}{%
\typeout{Selecting twosided page geometry}
\geometry{a4paper, left=30mm, right=30mm, top=20mm, bottom=35mm,
  footskip=20mm, footnotesep=\baselineskip, nohead}
}

% Load the titlesec package to control the formatting of titles
\RequirePackage[sf,bf,raggedright,noindentafter]{titlesec}[2002/04/07]
\ifthenelse{\boolean{longreport}}{%
  \titleformat{\part}[display]%
    {\centering\sffamily\bfseries\Huge}%
    {\partname\ \thepart}%
    {0pt}{}[\thispagestyle{empty}]
  \titleformat{\chapter}[hang]%
    {\sffamily\bfseries}%
    {\Huge\thechapter}%
    {1em}{\huge}[\thispagestyle{fancy}]
  \titlespacing{\chapter}{0pt}{-\baselineskip}{10ex}
  \titleformat{\section}[hang]%
    {\sffamily\bfseries}%
    {\Large\thesection}%
    {0.5em}{\Large}
  \titlespacing{\section}{0pt}{\baselineskip-\parskip}{-\parskip}
  \titleformat{\subsection}[hang]
    {\sffamily\bfseries}%
    {\normalsize\thesubsection}%
    {0.5em}{\normalsize}}{% else
  \titleformat{\section}[hang]%
    {\sffamily\bfseries}%
    {\LARGE\thesection}%
    {0.5em}{\LARGE}
  \titlespacing{\section}{0pt}{1.5\baselineskip-0.5\parskip}{\baselineskip-\parskip}
  \titleformat{\subsection}[hang]
    {\sffamily\bfseries}%
    {\large\thesubsection}%
    {0.5em}{\large}
  \titleformat{\subsubsection}[hang]
    {\sffamily\bfseries}%
    {\normalsize\thesubsubsection}%
    {0.5em}{\normalsize}}
\titlespacing{\subsection}{0pt}{\baselineskip-\parskip}{-\parskip}
\titlespacing{\subsubsection}{0pt}{\baselineskip-\parskip}{-\parskip}


\RequirePackage[absolute]{textpos}[2003/09/07]


% Format the captions
\RequirePackage{caption}[2004/01/23]
\captionsetup{%
  format=default,
  labelsep=period,
  singlelinecheck=false,
  font={small, sf}}

\RequirePackage{subfig}
\RequirePackage{eso-pic}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Miscellaneous formatting
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Turn off Victorian inter-sentence spacing.
\frenchspacing

% Fix layout of list environments with utilities from the paralist
% package.
\setdefaultitem{\textbullet}%\textperiodcentered
               {--}
               {\textasteriskcentered}%{\footnotesize\textbullet}
               {\textperiodcentered}
\setdefaultenum{1.}{a.}{1)}{a)}
\setdefaultleftmargin{2em}{2em}{2em}{2em}{1em}{1em}

% Fix footnote marker and rule
\def\@makefntext#1{\noindent\@makefootmark.\quad#1}
\def\@makefootmark{\@thefnmark}
\renewcommand{\footnoterule}{%
  \noindent\color[gray]{0.55}\rule{\textwidth}{0.5pt}\vspace{0.5\baselineskip}\normalcolor}



\providecommand{\clearemptydoublepage}{%
\ifthenelse{\boolean{screen}}{%
  \clearpage}{%
  \newpage{\pagestyle{empty}\cleardoublepage}}
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Front page definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\setlength{\TPHorizModule}{1mm}
\setlength{\TPVertModule}{1mm}
\textblockorigin{0mm}{0mm}

\providecommand{\frontpageprintelements}{%
  \ifthenelse{\boolean{report}}{%
    \ifthenelse{\boolean{norsk_selected} \or
      \boolean{nynorsk_selected}}{%
      \designelement{elements/rapport-s1}
    }{%
      \designelement{elements/report-s1}
    }}{%
    \ifthenelse{\boolean{norsk_selected} \or
      \boolean{nynorsk_selected}}{%
      \designelement{elements/notat-s1}
    }{%
      \designelement{elements/note-s1}
    }}
}

\providecommand{\frontpageprinttitle}{
    \begin{textblock}{122}(58,50)
      \noindent\begin{minipage}{122mm}
        \fontfamily{\sfdefault}
        \fontseries{bx}
        \fontsize{30}{36}
        \selectfont
        \raggedright
        \noindent\@title
        \ifthenelse{\equal{\@subtitle}{}}{}{
          \vspace{0.5em} 
          \fontsize{20}{24} \selectfont\newline 
          \@subtitle
        }        
      \end{minipage}
    \end{textblock}
}
\providecommand{\frontpageprintmatter}{%
  \begin{textblock}{185}(18,234)
    \parindent=\z@
    \begin{sffamily}
      \bfseries
      \begin{tabular}{@{}>{\raggedleft}p{25.5mm}@{}p{14.5mm}@{}>{\raggedright\let\\\tabularnewline}p{122mm}}
        \ifthenelse{\boolean{report}}{%
          \colorbox{nrblue}{\textcolor{white}{\reportnumbername}}
        }{%
          \colorbox{nrblue}{\textcolor{white}{\notenumbername}}
        }
        && \printreportnumber\\
        \colorbox{nrblue}{\textcolor{white}{\authorname}}
        && 
        \ifthenelse{\value{numberofauthors}>\value{largenumberofauthors}}{\vspace{-2em}\begin{multicols}{2}
        }{}
        \@author 
        \ifthenelse{\value{numberofauthors}>\value{largenumberofauthors}}{\end{multicols}}{}
        \\
        \colorbox{nrblue}{\textcolor{white}{\datename}}
        && \@date\\
        \ifthenelse{\boolean{report}}{%
          \colorbox{nrblue}{\textcolor{white}{\isbnname}}}{}
        && \ifthenelse{\boolean{report}}{\printisbn}{}
      \end{tabular}
    \end{sffamily}
  \end{textblock}
}
\providecommand{\frontpageprintfigure}{%
  \begin{textblock}{147}(58,114)
    \noindent
    \ifthenelse{\boolean{frontpagefiguredefined}}{%
      \includegraphics[width=120mm]{\frontpagefigurefile}%
      \ifthenelse{\equal{\frontpagefigurecpright}{}}{}{\newline \small \copyright\ \frontpagefigurecpright}}%
     {\frontpagefigurecommands}
  \end{textblock}
}

\providecommand{\nrfrontpage}{%
  \ifthenelse{\boolean{draftversion}}{
    \AddToShipoutPicture*{\frontpageprintelements}	
    \thispagestyle{empty}%
    \noindent
    \setboolean{hidden}{true}%\typeout{Making front page!}
        % Write authors outside page.    
    \makebox[0pt]{\mbox{}\hspace{100cm}\textcolor{white}{\@author}}
    \setboolean{hidden}{false}
    \frontpageprinttitle
    \frontpageprintmatter
    \frontpageprintfigure
    % Trick LaTeX into believing there is something here...
    \hspace{0pt}
    \clearpage
   }
   {% else final version
    \AddToShipoutPicture*{\frontpageprintelements}
    \thispagestyle{empty}%
    \noindent
    \setboolean{hidden}{true}%\typeout{Making front page!}
        % Write authors outside page.    
    \makebox[0pt]{\mbox{}\hspace{100cm}\textcolor{white}{\@author}}
    \setboolean{hidden}{false}
    \frontpageprinttitle
    \frontpageprintmatter
    \frontpageprintfigure
    % Trick LaTeX into believing there is something here...
    \hspace{0pt}
    \clearpage
  }
}


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   ``About'' page definition (page 2)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\providecommand{\aboutpageentry}[2]{%
  \noindent\textbf{\sffamily\large#1}\\
  #2}

\providecommand{\aboutauthorstext}{%
Use \texttt{$\backslash$aboutauthors} to enter author information.}

\providecommand{\aboutnrtext}{%
Norsk Regnesentral (Norwegian Computing Center, NR) is a private,
independent, non-profit foundation established in 1952. NR carries out
contract research and development projects in information and
communication technology and applied statistical-mathe\-matical
modelling. The clients include a broad range of industrial, commercial
and public service organisations in the national as well as the
international market. Our scientific and technical capabilities are
further developed in co-operation with The Research Council of Norway
and key customers. The results of our projects may take the form of
reports, software, prototypes, and short courses. A proof of the
confidence and appreciation our clients have in us is given by the
fact that most of our new contracts are signed with previous
customers. }

% A separate translation into nynorsk should be provided here
\ifthenelse{\boolean{norsk_selected} \or
  \boolean{nynorsk_selected}}{%
  \renewcommand{\aboutnrtext}{%
Norsk Regnesentral (NR) er en privat, uavhengig stiftelse som utf{\o}rer
oppdrags\-forskning for bedrifter og det offentlige i det norske og
internasjonale markedet. NR ble etablert i 1952 og har kontorer i
Kristen Nygaards hus ved Universitetet i Oslo. NR er et av Europas
st{\o}rste milj{\o}er innen anvendt statistisk-matematisk modellering og har
et senter for forskningsdrevet innovasjon, Statistics for Innovation --
(sfi)$^2$, med finansiering fra Norges forskningsr{\aa}d. Det jobbes med et
bredt spekter av problemstillinger, for eksempel finansiell risiko,
jord\-observasjon, estimering av fiske\-bestander og beskrivelse av
geologien i petroleums\-reservoarer. NR er ledende i Norge innen
utvalgte deler av informasjons- og kommunikasjonsteknologi. Innen
IKT-omr{\aa}det har NR innsatsomr{\aa}dene e-inkludering,
infor\-ma\-sjons\-sikkerhet og smarte infor\-masjons\-systemer. 

\begin{center}
NRs visjon er forskningsresultater som brukes og synes.
\end{center}
}}{}
% The command \aboutothers is provided to allow entries for other
% organisations, eg. if an author is external to NR. A
% redefinition should be of the form:
%  \renewcommand{\aboutothers}{%
%    \aboutpageentry{Company A}{%
%      Company A specialises in ...}}
% If several entries are required, use multiple calls to
% \aboutpageentry and separate them with \par\bigskip as in
% \aboutpagetext below.
\providecommand{\aboutotherstext}{}

\providecommand{\aboutothers}[2]{%
\renewcommand{\aboutotherstext}{\aboutpageentry{#1}{#2}}}

\providecommand{\printaboutpagetext}{%
  \ifthenelse{\boolean{aboutauthors}}{%
  \aboutpageentry{\aboutauthorsheading}{\aboutauthorstext}\par\bigskip}{}
  \aboutpageentry{\nrname}{\aboutnrtext}\par\bigskip
  \aboutotherstext\par\bigskip
}


\providecommand{\nraboutpage}{%
  \thispagestyle{empty}
  \noindent
  \designelement{elements/NR-logo_utvidet_60-5mm_cmyk}
  \par\vspace{4\baselineskip}
  \printaboutpagetext
  \par\vfill\noindent
  \begin{textblock}{151}(30,267)
    \noindent
    \designelement{elements/NR-address}
  \end{textblock}
  \clearpage
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Title page definition
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\newdimen\@bls   \@bls\baselineskip
\newdimen\@ls    \@ls\lineskip
\newdimen\@lsl   \@lsl\lineskiplimit
\def\@oninterlineskip{\baselineskip\@bls\lineskip\@ls\lineskiplimit\@lsl}
\newenvironment{titlelist}{\offinterlineskip\begin{list}{}%
    {\labelwidth=1.5em%
     \labelsep=0pt%
     \leftmargin=0pt%
     \rightmargin=0pt%
     \topsep=0pt%
     \parsep=0pt%
     \partopsep=0pt%
     \itemsep=0.5\@bls%
     \sffamily}}{\end{list}}
\def\t@textbox#1#2{\begin{minipage}[t]{#1}\@oninterlineskip\raggedright
    \strut#2\strut
  \end{minipage}}
\def\t@text#1#2{\t@textbox{0.3\textwidth}{#1}\hfill\t@textbox{0.67\textwidth}{#2}}
\def\titleitem#1#2{\item\t@text{#1}{#2}}

\providecommand{\nrtitlepage}{%
  \thispagestyle{plain}
  \setboolean{titlepage}{true}
  \renewcommand{\and}{,\ }
  %\parindent\z@
  \begin{titlelist}
    \titleitem{\textbf{\large\titlename}}{\textbf{\large\printtitle}}
    \titleitem{\textbf{\authorname}}     {\textbf{\@author}}
    \printqatext
   \titleitem{\datename}                {\@date}
%    \titleitem{\yearname}                {\printyear}
    \ifthenelse{\boolean{report}}{\titleitem{\isbnname}{\printisbn}}{}
    \titleitem{\publicationnumbername}   {\printreportnumber}
  \end{titlelist}
  \par\vspace{4\baselineskip}
}

\newboolean{titlepagemade}

% Mainmatter is not needed (kept for backwards compatability)
\ifthenelse{\boolean{longreport}}{
  \renewcommand{\mainmatter}{\typeout{NRdoc info: The command mainmatter is not needed.}}}{%
  \providecommand{\mainmatter}{\typeout{NRdoc info: The command mainmatter is not needed.}}}

\ifthenelse{\boolean{longreport}}{%
  \newenvironment{abstract}{%
    \ifthenelse{\not\boolean{titlepagemade}}{%
      \typeout{NRdoc warning: abstract must be put after "\protect\titlepage" command.}}{}
    \noindent\textbf{\sffamily\large\abstractname}\\*}
  {\vfill
    %\parindent\z@
    \begin{titlelist}
      \titleitem{\keywordsname}        {\printkeywords}
      \titleitem{\targetgroupname}     {\printtarget}
      \titleitem{\availabilityname}    {\printavailability}
      \titleitem{\projectname}         {\printproject}
      \titleitem{\projectnumbername}   {\printprojectnumber}
      \titleitem{\researchfieldname}   {\printresearchfield}
      \titleitem{\numberofpagesname}   {\pageref*{LastPage}}  

      \titleitem{Copyright \copyright\ \printyear}{\nrname}
    \end{titlelist}
  }}{%
  \renewenvironment{abstract}{%
    \ifthenelse{\not\boolean{titlepagemade}}{%
      \typeout{NRdoc warning: abstract must be put after "\protect\titlepage" command.}}{}
    \noindent\textbf{\sffamily\large\abstractname}\\*}
  {\vfill\noindent
    %\parindent\z@
    \begin{titlelist}
      \titleitem{\keywordsname}        {\printkeywords}
      \titleitem{\targetgroupname}     {\printtarget}
      \titleitem{\availabilityname}    {\printavailability}
      \titleitem{\projectname}         {\printproject}
      \titleitem{\projectnumbername}   {\printprojectnumber}
      \titleitem{\researchfieldname}   {\printresearchfield}
      \titleitem{\numberofpagesname}   {\pageref*{LastPage}}    
      \titleitem{\copyright\ Copyright}{\nrname}
    \end{titlelist}
  }
  \renewcommand\appendix{\par
  \clearpage%
  \setcounter{section}{0}%
  \setcounter{subsection}{0}%  
  \gdef\thesection{\@Alph\c@section}%
  \gdef\theHsection{\@Alph\c@section}%
  \setcounter{table}{0}% 
  \renewcommand{\thetable}{\thesection.\arabic{table}}%  
  \gdef\theHtable{\theHsection.\arabic{table}}%
  \setcounter{equation}{0}%
  \renewcommand{\theequation}{\thesection.\arabic{equation}}%
  \gdef\theHequation{\theHsection.\arabic{equation}}%
  \setcounter{figure}{0}%
  \renewcommand{\thefigure}{\thesection.\arabic{figure}}%
  \gdef\theHfigure{\theHsection.\arabic{figure}}%
  \gdef\theHsubfigure{\theHfigure.\arabic{subfigure}}%
  }%
}

\renewcommand{\maketitle}{%
  \setboolean{titlepagemade}{true}
  \hypersetup{pdftitle={\printtitle},
              pdfkeywords={\pdfprintkeywords},
              pdfsubject={\printsubject},
              pdfauthor={\printauthorforpdf}
        }% 
  %\setboolean{pdfdataset}{true}
  \pagenumbering{arabic}
  \nrfrontpage%  
  \nraboutpage
  \nrtitlepage  
  \relax  
}




%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Toc, lof etc.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\renewcommand{\tableofcontents}{%
  \clearemptydoublepage
  \ifthenelse{\boolean{longreport}}{%
    \chapter*{\contentsname}}{%
    \section*{\contentsname}}
  \addtocontents{toc}{\protect\hypertarget{toc}{}%
    \protect\pdfbookmark[1]{\contentsname}{toc}}
  \@starttoc{toc}% 
  \clearemptydoublepage
}

\renewcommand{\listoffigures}{%
  \clearemptydoublepage
  \ifthenelse{\boolean{longreport}}{%
    \chapter*{\listfigurename}}{%
    \section*{\listfigurename}}
  \addtocontents{lof}{\protect\hypertarget{lof}{}%
    \protect\pdfbookmark[1]{\listfigurename}{lof}}
  \@starttoc{lof}%
  \clearemptydoublepage
}

\renewcommand{\listoftables}{% 
  \clearemptydoublepage
  \ifthenelse{\boolean{longreport}}{%
    \chapter*{\listtablename}}{%
    \section*{\listtablename}}
  \addtocontents{lot}{\protect\hypertarget{lot}{}%
    \protect\pdfbookmark[1]{\listtablename}{lot}}
  \@starttoc{lot}%
  \clearemptydoublepage
}

\newcommand{\listofexamples}{% 
  \clearemptydoublepage
  \ifthenelse{\boolean{longreport}}{%
    \chapter*{\listexamplename}}{%
    \section*{\listexamplename}}
  \addtocontents{loe}{\protect\hypertarget{loe}{}%
    \protect\pdfbookmark[1]{\listexamplename}{loe}}
  \@starttoc{loe}%
  \clearemptydoublepage
}


% Example environment
\ifthenelse{\boolean{longreport}}{
  \newcounter{example}[chapter]
  \def\theexample{\thechapter.\@arabic\c@example}
}{
  \newcounter{example}[section]
  \def\theexample{\thesection.\@arabic\c@example}
}
\let\l@example\l@figure
\providecommand\ext@example{loe}
\providecommand\fnum@example{\examplename~\theexample}
\def\toclevel@example{0}
\def\ftype@example{2}

\newenvironment{example}{%

  \list{}{\topsep 0pt%
          \itemsep 0pt%
          \labelwidth 0pt%
          \leftmargin 0pt%
          \listparindent \parindent%
          \itemindent 0pt
          \parsep 0pt}\item\relax%  
 \def\@captype{example}%
}{\endlist}

\captionsetup[example]{%
  belowskip=10pt,
  position=above,
  format=default,
  labelsep=period,  
  singlelinecheck=false,
  font={normalsize, rm},
  labelfont=bf%
}

% Exercise
\newcommand{\listofexercises}{% 
  \clearemptydoublepage
  \ifthenelse{\boolean{longreport}}{%
    \chapter*{\listexercisename}}{%
    \section*{\listexercisename}}
  \addtocontents{loe}{\protect\hypertarget{loe}{}%
    \protect\pdfbookmark[1]{\listexercisename}{loe}}
  \@starttoc{loe}%
  \clearemptydoublepage
}


% Exercise environment
\ifthenelse{\boolean{longreport}}{
  \newcounter{exercise}[chapter]
  \def\theexercise{\thechapter.\@arabic\c@exercise}
}{
  \newcounter{exercise}[section]
  \def\theexercise{\thesection.\@arabic\c@exercise}
}
\let\l@exercise\l@figure
\providecommand\ext@exercise{loe}
\providecommand\fnum@exercise{\exercisename~\theexercise}
\def\toclevel@exercise{0}
\def\ftype@exercise{2}

\newenvironment{exercise}{%

  \list{}{\topsep 0pt%
          \itemsep 0pt%
          \labelwidth 0pt%
          \leftmargin 0pt%
          \listparindent \parindent%
          \itemindent 0pt
          \parsep 0pt}\item\relax%  
 \def\@captype{exercise}%
}{\endlist}

\captionsetup[exercise]{%
  belowskip=10pt,
  position=above,
  format=default,
  labelsep=period,  
  singlelinecheck=false,
  font={normalsize, rm},
  labelfont=bf%
}

% Titletoc formats the table of contents, LoF etc.
\RequirePackage{titletoc}[2002/03/27]
\contentsmargin{2em}
\ifthenelse{\boolean{longreport}}{
  \titlecontents{part}[2em]%
    {\bigskip\sffamily\bfseries}{\contentslabel{2em}}%
    {\hspace{-2em}}%
    {\titlerule*[1pc]{}\contentspage}
  \titlecontents{chapter}[2em]%
    {\bigskip\sffamily\bfseries}{\contentslabel{2em}}%
    {\hspace{-2em}}%
    {\textnormal{\sffamily\titlerule*[1pc]{.}}\contentspage}
  \titlecontents{section}[5em]%
    {\sffamily}{\contentslabel{3em}}%
    {\hspace{-3em}}%
    {\titlerule*[1pc]{.}\contentspage}
  \titlecontents{subsection}[9em]%
    {\sffamily}{\contentslabel{4em}}%
    {\hspace{-4em}}%
    {\titlerule*[1pc]{.}\contentspage}}{% else
  \titlecontents{section}[2em]%
    {\bigskip\sffamily\bfseries}{\contentslabel{2em}}%
    {\hspace{-2em}}%
    {\textnormal{\sffamily\titlerule*[1pc]{.}}\contentspage}
  \titlecontents{subsection}[5em]%
    {\sffamily}{\contentslabel{3em}}%
    {\hspace{-3em}}%
    {\titlerule*[1pc]{.}\contentspage}
  \titlecontents{subsubsection}[9em]%
    {\sffamily}{\contentslabel{4em}}%
    {\hspace{-4em}}%
    {\titlerule*[1pc]{.}\contentspage}}
\titlecontents{figure}[2em]%
  {\sffamily}{\contentslabel{2em}}%
  {\hspace{-2em}}%
  {\textnormal{\sffamily\titlerule*[1pc]{.}}\contentspage}
\titlecontents{table}[2em]%
  {\sffamily}{\contentslabel{2em}}%
  {\hspace{-2em}}%
  {\textnormal{\sffamily\titlerule*[1pc]{.}}\contentspage}
\contentsuse{example}{loe}
\titlecontents{example}[2em]%
  {\sffamily}{\contentslabel{2em}}%
  {\hspace{-2em}}%
  {\textnormal{\sffamily\titlerule*[1pc]{.}}\contentspage}
\contentsuse{exercise}{loe}
\titlecontents{exercise}[2em]%
  {\sffamily}{\contentslabel{2em}}%
  {\hspace{-2em}}%
  {\textnormal{\sffamily\titlerule*[1pc]{.}}\contentspage}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Fancy footers
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\pagestyle{fancy}
\renewcommand{\headrulewidth}{0pt}
\renewcommand{\footrulewidth}{0pt}
\newsavebox{\footerlogo}
\savebox{\footerlogo}{\raisebox{-0.5ex}{\designelement{elements/NR-logo_footer}}}
\ifthenelse{\boolean{@twoside}}
{%
  \fancyhead{}
  \fancyfoot{}
  \fancyfoot[RO]{\sffamily\bfseries\footnotesize%
    \printshorttitle{}%
    \qquad\usebox{\footerlogo}%
    \qquad\thepage}
  \fancyfoot[LE]{\sffamily\bfseries\footnotesize%
    \thepage%
    \qquad\usebox{\footerlogo}%
    \qquad\printshorttitle{}}
}{%
  \fancyhead{}
  \fancyfoot{}
  \fancyfoot[RO]{\sffamily\bfseries\footnotesize%
    \printshorttitle{}%
    \qquad\usebox{\footerlogo}%
    \qquad\thepage}
}
% Fancy footer for ``plain'' pages also
\ifthenelse{\boolean{@twoside}}
{%
\fancypagestyle{plain}{%
  \fancyhead{}
  \fancyfoot{}
  \fancyfoot[RO]{\sffamily\bfseries\footnotesize%
    \qquad\usebox{\footerlogo}%
    \qquad\thepage}
  \fancyfoot[LE]{\sffamily\bfseries\footnotesize%
    \thepage%
    \qquad\usebox{\footerlogo}}
}
}{%
\fancypagestyle{plain}{%
  \fancyhead{}
  \fancyfoot{}
  \fancyfoot[RO]{\sffamily\bfseries\footnotesize%
    \qquad\usebox{\footerlogo}%
    \qquad\thepage}
 
}       
}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Bibliography and index
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\ifthenelse{\boolean{nobibstyle}}{%
  \typeout{NRdoc warning: No bibliography style set  by package. Use your own}
}%
{
  \ifthenelse{\boolean{citenumeric}}{%
    \ifthenelse{\boolean{citealphanumeric}}{%
      \ifthenelse{\boolean{norsk_selected}%
        \or%
        \boolean{nynorsk_selected}}{%    
        \bibliographystyle{apalike-url-norsk}
      }{%        
        \bibliographystyle{apalike-url}
      }
    }{%
      \bibliographystyle{unsrturl} %\bibliographystyle{nrunsrt}
    }
    \typeout{nrdoc uses numeric citation labels}
  }{%
    \ifthenelse{\boolean{norsk_selected}%
      \or%
      \boolean{nynorsk_selected}}{%
      \typeout{nrdoc uses alphabethic citation labels}
      \bibliographystyle{apalike-url-norsk}
    }{%
      \bibliographystyle{apalike-url} %\bibliographystyle{nrplain}
      \typeout{nrdoc uses alphabethic citation labels}
    }
  }
}

\renewenvironment{thebibliography}[1]{% 
  \ifthenelse{\boolean{longreport}}{%
    \phantomsection
    \chapter*{\refname}\addcontentsline{toc}{chapter}{\refname}}{%
    \phantomsection
    \section*{\refname}\addcontentsline{toc}{section}{\refname}}
  \list{\@biblabel{\@arabic\c@enumiv}}{%
    \settowidth\labelwidth{\@biblabel{#1}}%
    \leftmargin\labelwidth
    \advance\leftmargin\labelsep
    \@openbib@code
    \usecounter{enumiv}%
    \let\p@enumiv\@empty
    \renewcommand\theenumiv{\@arabic\c@enumiv}}%
  \sloppy
  \clubpenalty4000
  \@clubpenalty \clubpenalty
  \widowpenalty4000%
  \sfcode`\.\@m}%
{\def\@noitemerr
  {\@latex@warning{Empty `thebibliography' environment}}%
  \endlist}

\renewenvironment{theindex}{%
  \clearemptydoublepage
  \ifthenelse{\boolean{longreport}}{% 
    \phantomsection
    \chapter*{\indexname}%
    \addcontentsline{toc}{chapter}{\indexname}}{%
    \phantomsection
    \section*{\indexname}%
    \addcontentsline{toc}{section}{\indexname}}
  \begin{multicols}{2}
  \thispagestyle{fancy}%\parindent\z@
  \parskip\z@ \@plus .3\p@\relax
  \let\item\@idxitem}
  {\end{multicols}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   LastPage
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\def\lastpage@putlabel{\addtocounter{page}{-1}%
  \immediate\write\@auxout{\string
    \newlabel{LastPage}{{}{\thepage}{\relax}{}{}}}%
  \addtocounter{page}{1}}

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%   Misc
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

\ifthenelse{\boolean{longreport}}{%
  \renewcommand\chapter{%
    \clearemptydoublepage
    \global\@topnum\z@
    \@afterindentfalse
    \secdef\@chapter\@schapter}}{}

\providecommand{\floor}[1]{\lfloor #1 \rfloor}
\providecommand{\ceil}[1]{\lceil #1 \rceil}

\ifpdf
\DeclareGraphicsExtensions{%
.jpg,.JPG,.jpeg,.JPEG,.png,.PNG,.tif,.TIF,.tiff,.pdf,.PDF}
\DeclareGraphicsRule{.JPG}{jpg}{.JPG}{}
\DeclareGraphicsRule{.JPEG}{jpg}{.JPEG}{}
\DeclareGraphicsRule{.jpeg}{jpg}{.jpeg}{}
\DeclareGraphicsRule{.PNG}{png}{.PNG}{}
\DeclareGraphicsRule{.TIF}{tif}{.TIF}{}
\DeclareGraphicsRule{.tiff}{tif}{.tiff}{}
\DeclareGraphicsRule{.PDF}{pdf}{.PDF}{}
\else 
\DeclareGraphicsExtensions{.eps,.EPS,.ps,.PS}
\DeclareGraphicsRule{.EPS}{eps}{.EPS}{}
\DeclareGraphicsRule{.PS}{ps}{.PS}{}
\fi


\AtBeginDocument{\renewcommand{\l@table}{\@dottedtocline{1}{0em}{3.5em}}
\renewcommand{\l@figure}{\@dottedtocline{1}{0em}{3.5em}} 
}

\AtEndDocument{%
  \ifthenelse{\boolean{noindex}}%
        {%
        }{%
        \printindex%
        }
  \clearpage\lastpage@putlabel%
}
