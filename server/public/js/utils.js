// Page utils
export class PageUtils {
    static clearFormFields() {
        $('form').find('input[type=text], input[type=number], input[type=email], input[type=password], select').val('');
    }

    static clickOnEnter(targetBtnId) {
        document.addEventListener('keydown', function(event) {
            if (event.key === 'Enter') {
                event.preventDefault();
                document.getElementById(targetBtnId).click();
            }
        });
    }

    static startSpinner(button) {
        button.innerHTML = '<i class="fas fa-spinner fa-spin"></i>';
        button.setAttribute('disabled', true);
    }

    static stopSpinner(button,text) {
        button.innerHTML = text;
        button.removeAttribute('disabled');
    }

    static toggleVisibility(emitterSelector, receiverSelector, condition) {
        $(emitterSelector).click(function() {
            const emitter = this; //DOM object
            const receiver = $(receiverSelector).get(0); //DOM object
    
            const shouldShow = condition(emitter, receiver);
    
            $(receiver).toggleClass('d-none', !shouldShow);
            $(receiver).val('');
            $(receiver).find('input, select').val('');
        });
    }  

    static showToast(header, message, color = 'yellow') {
        // Check if the toast container exists
        if ($('#submissionToast').length === 0) {
            // Toast HTML structure
            const toastHTML = `
                <div style="position: fixed; top: 80px; right: 20px; min-height: 200px">
                    <div style="position: absolute; top: 0; right: 0;">
                        <div id="submissionToast" class="toast" role="alert">
                            <div class="toast-header">
                                <strong class="mr-auto" style="color: black">Submission Status</strong>
                                <button type="button" class="ml-2 mb-1 close" data-dismiss="toast">
                                    <span >&times;</span>
                                </button>
                            </div>
                            <div class="toast-body" style="width: 400px;"></div>
                        </div>
                    </div>
                </div>
            `;

            // Append the toast HTML to the body
            $('body').append(toastHTML);

            // Initialize the toast (if you're using a library that requires initialization)
            $('#submissionToast').toast({ delay: 3000 }); // Example for Bootstrap toast
        }

        // Set colors based on the 'color' parameter
        const colors = color === 'red' ? 
            {header: 'lightcoral', body: 'lightpink'} : 
            color === 'green' ? 
                {header: 'rgb(173, 234, 173)', body: 'lightcyan'} :
                {header: 'khaki', body: 'khaki'};

        // Update the toast content and colors
        $('#submissionToast .toast-header strong').html(header);
        $('#submissionToast .toast-body').text(message).css('background-color', colors.body);
        $('#submissionToast').find('.toast-header').css('background-color', colors.header);

        // Show the toast
        $('#submissionToast').toast('show');
    }

    static downloadTxtFile(text) {
        const element = document.createElement('a');
        const file = new Blob([text], {type: 'text/plain'});
        element.href = URL.createObjectURL(file);
        element.download = "download_" + new Date().toISOString().replace(/[\W_]+/g, "_") + ".txt";
        document.body.appendChild(element); // Required for this to work in FireFox
        element.click();
        document.body.removeChild(element);
    }
};

export class UserUtils { 
    static readCookie(name) {
        const nameEQ = name + "=";
        const ca = document.cookie.split(';');
        for(let i=0;i < ca.length;i++) {
            let c = ca[i];
            while (c.charAt(0)==' ') c = c.substring(1,c.length);
            if (c.indexOf(nameEQ) == 0) return c.substring(nameEQ.length,c.length);
        }
        return null;
    }

    static displayUserAlias(divId) {
        const alias = UserUtils.getUserAlias();
        if (alias === 'Anonymous') {
            document.getElementById(divId).innerHTML = '<a class="btn btn-success" href="/login">Sign In</a>';
        } else {
            document.getElementById(divId).innerHTML = `<i class="fa-solid fa-circle-user"></i> <strong>${alias}</strong>`;
        }
    };

    static getUserAlias() {
        const userDataCookie = UserUtils.readCookie('userData');
        if (userDataCookie) {
            const userData = JSON.parse(decodeURIComponent(userDataCookie));
            return userData.alias;
        } else {
            return 'Anonymous';
        }
    }

    static getUserId() {
        const userDataCookie = UserUtils.readCookie('userData');
        if (userDataCookie) {
            const userData = JSON.parse(decodeURIComponent(userDataCookie));
            return userData._id;
        } else {
            return null;
        }
    }
};

export class HTTPUtils {
    static async fetchUniprotInfo(accessionNumbers) {
        const accessionNumberList = accessionNumbers.split(',').map(an => an.trim());
    
        let uniprotPdbCodes = new Set();
        let mass = [];
        let length = [];
    
        for (const accessionNumber of accessionNumberList) {
            const url = `https://rest.uniprot.org/uniprotkb/${accessionNumber}.json?fields=xref_pdb,length,mass`;
    
            try {
                const response = await fetch(url);
                if (!response.ok) {
                    throw new Error(`HTTP error! status: ${response.status}`);
                }
                const data = await response.json();
    
                if (data.uniProtKBCrossReferences) {
                    data.uniProtKBCrossReferences
                        .filter(xref => xref.database === "PDB")
                        .forEach(xref => uniprotPdbCodes.add(xref.id));
                }
    
                mass.push(data.sequence ? data.sequence.molWeight : null);
                length.push(data.sequence ? data.sequence.length : null);
                
            } catch (error) {
                console.error('Error fetching data for accession number:', accessionNumber, error);
            }
        }
    
        uniprotPdbCodes = Array.from(uniprotPdbCodes);
        return { uniprotPdbCodes, mass, length };
    }
}

